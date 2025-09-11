# Frontend Application Structure and Architecture

## Overview
This front-end React/Vite dashboard is a comprehensive scientific application organized into modular components for ubiquitin chain visualization, synthesis planning, and structural analysis. The architecture follows a hierarchical structure with clear separation between application-level configuration, styling, and component-specific functionality.

---

## Application Foundation

### 1. `main.jsx` - Application Entry Point
**Role:** React application bootstrapper and root renderer.

**Technical Implementation:**
```jsx
import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.jsx'

createRoot(document.getElementById('root')).render(
  <StrictMode>
    <App />
  </StrictMode>,
)
```

**Key Features:**
- **StrictMode Wrapper:** Enables React's development-time checks and warnings
- **Root Element Mounting:** Connects React application to DOM element with id "root"
- **CSS Import:** Loads global styling before application initialization
- **Modern React API:** Uses React 18's createRoot API for improved performance

### 2. `App.jsx` - Application Shell
**Role:** Top-level application component and layout manager.

**Structure & Styling:**
```jsx
function App() {
  return (
    <div className="App">
      <h1 style={{ 
        textAlign: 'center', 
        margin: '32px 0 24px 0', 
        color: 'black', 
        fontWeight: 'bold', 
        textShadow: '0 2px 8px #fff' 
      }}>
        Ubiquitin Multimer Synthesis Dashboard
      </h1>
      <ModuleDashboard />
    </div>
  );
}
```

**Design Decisions:**
- **Minimal Shell:** Provides application title and container for ModuleDashboard
- **Inline Styling:** Strategic use of inline styles for the main heading
- **Accessibility:** Clear hierarchy with semantic h1 element
- **Single Component Import:** Delegates all functionality to ModuleDashboard

### 3. `index.css` - Global Styling Foundation
**Role:** Establishes global styling foundation and design system variables.

**Core Design System:**
```css
:root {
  font-family: system-ui, Avenir, Helvetica, Arial, sans-serif;
  line-height: 1.5;
  font-weight: 400;
  color-scheme: light dark;
  color: rgba(255, 255, 255, 0.87);
  background-color: #242424;
}
```

**Key Features:**
- **Dark Theme Default:** Primary background #242424 with light text
- **System Font Stack:** Modern, cross-platform font selection
- **Color Scheme Support:** Automatic light/dark mode detection
- **Responsive Typography:** Optimized line-height and font-weight
- **Cross-Browser Optimization:** Font smoothing and rendering optimizations

**Interactive Element Styling:**
```css
button {
  border-radius: 8px;
  border: 1px solid transparent;
  padding: 0.6em 1.2em;
  font-size: 1em;
  font-weight: 500;
  font-family: inherit;
  background-color: #1a1a1a;
  cursor: pointer;
  transition: border-color 0.25s;
}
button:hover {
  border-color: #646cff;
}
```

**Design Principles:**
- **Consistent Interactivity:** Smooth transitions and hover effects
- **Accessibility Focus:** Focus-visible states for keyboard navigation
- **Theme Flexibility:** Light mode overrides with `@media (prefers-color-scheme: light)`
- **Performance Optimized:** Hardware-accelerated animations

### 4. `App.css` - Component-Specific Styling
**Role:** Provides component-specific styling and animation utilities.

**Logo Animation System:**
```css
.logo {
  height: 6em;
  padding: 1.5em;
  will-change: filter;
  transition: filter 300ms;
}

@keyframes logo-spin {
  from { transform: rotate(0deg); }
  to { transform: rotate(360deg); }
}

@media (prefers-reduced-motion: no-preference) {
  a:nth-of-type(2) .logo {
    animation: logo-spin infinite 20s linear;
  }
}
```

**Accessibility Features:**
- **Motion Preferences:** Respects user's reduced motion preferences
- **Performance Optimization:** Uses `will-change` for animation performance
- **Smooth Transitions:** 300ms filter transitions for visual feedback
- **Responsive Design:** Relative sizing with em units

**Layout Utilities:**
- **Max-width Container:** Centered layout with 1280px maximum width
- **Responsive Padding:** 2rem padding for comfortable spacing
- **Center Alignment:** Text-align center for primary content
- **Card Styling:** Reusable .card class with 2em padding

---

## Component Architecture

### 5. `ModuleDashboard.jsx` - Main Dashboard Controller
**Role:** Central application controller and page router.

**Key Features:**
- **Page Management:** Manages 6 pages through PAGE_CONFIG: draw, tetramers, pentamers, reactionStats, subgraph, nomenclature
- **State Management:** Coordinates state across all child components including:
  - Selected panels for synthesis planning
  - API response data (figures, reaction sequences, nomenclature data)
  - Interactive scaffold state (nodes, arrows, mapping)
- **API Integration:** Handles all backend communication for synthesis planning and nomenclature analysis
- **Dynamic Rendering:** Conditionally renders different page layouts based on current page selection

**Key Imports & Dependencies:**
- `Panel` - Container styling
- `JsonToScaffold` - JSON-to-visual conversion
- `ClickableScaffoldPanel` - Interactive drawing interface
- `ReactionSequencesPaneled` - Synthesis pathway visualization
- `SubgraphAnalysisPage` - Isomorphism analysis
- `ReactionPathStatisticsPage` - Metrics visualization
- `EdgeTreeViewer` - Nomenclature tree display

---

## Page-Specific Components

### 6. `ClickableScaffoldPanel.jsx` - Interactive Drawing Interface
**Role:** Primary interactive component for the "Explore Reaction Pathways" page.

**Core Functionality:**
- **Canvas-Based Drawing:** Renders interactive node-edge graphs on HTML5 canvas
- **Click Logic:** Implements complex click-to-arrow logic for pathway construction
- **JSON Simulation:** Includes `simulateClicksFromJson()` function that converts JSON ubiquitin structures to visual representations
- **State Synchronization:** Maintains nodes, arrows, clickedNodes, and smallNodes state
- **API Integration:** Submits constructed pathways to backend for analysis

**Technical Details:**
- Uses refs for performance optimization in canvas operations
- Supports both interactive and frozen modes
- Handles complex graph traversal algorithms for pathway visualization
- Exports `simulateClicksFromJson` for use by other components

### 7. `JsonToScaffold.jsx` - Static Scaffold Renderer
**Role:** Converts JSON ubiquitin data into static visual representations.

**Key Features:**
- **JSON Processing:** Implements its own `simulateClicksFromJson()` function (similar to ClickableScaffoldPanel)
- **Preorder Traversal:** Follows biological conventions for ubiquitin chain numbering
- **Visual Logic:** K63 linkages move rightward/up, K48 linkages move leftward/up
- **Protection Groups:** Renders SMAC/ABOC protecting groups as small colored nodes
- **Static Display:** Uses ScaffoldPanel in frozen mode for non-interactive display

**Usage Context:**
- Used in tetramer/pentamer selection grids
- Renders individual multimer thumbnails
- Provides consistent visual representation across the application

### 8. `ScaffoldPanel.jsx` - Core Visualization Engine
**Role:** Low-level canvas rendering engine used by both interactive and static components.

**Technical Implementation:**
- **Canvas Management:** Handles all drawing operations (nodes, edges, arrows, small nodes)
- **State Management:** Tracks node states, arrow configurations, and small node positions
- **Event Handling:** Manages mouse interactions when not frozen
- **Visual Consistency:** Maintains consistent styling and colors across all scaffold displays
- **Modular Design:** Supports both interactive (ClickableScaffoldPanel) and static (JsonToScaffold) use cases

### 9. `ReactionSequencesPaneled.jsx` - Synthesis Pathway Display
**Role:** Visualizes step-by-step reaction sequences for synthesis planning.

**Features:**
- **Multi-Step Visualization:** Shows progression from dimer → trimer → tetramer → pentamer
- **Reaction Well Display:** Optional display of specific reaction conditions
- **JSON Processing:** Handles complex nested JSON structures from backend
- **Dynamic Layout:** Adapts to both tetramer and pentamer synthesis pathways
- **Integration:** Uses JsonToScaffold for rendering individual synthesis steps

### 10. `ReactionPathStatisticsPage.jsx` - Metrics Dashboard
**Role:** Statistical analysis and visualization for reaction pathways.

**Capabilities:**
- **Chart Integration:** Uses Chart.js for bar charts and scatter plots
- **API Communication:** Fetches statistical data from backend endpoints
- **Interactive Controls:** Allows filtering by multimer size (4 or 5) and pathway type (ABOC/SMAC)
- **Data Processing:** Handles complex statistical datasets and presents them visually

### 11. `SubgraphAnalysisPage.jsx` - Isomorphism Analysis
**Role:** Advanced analysis of structural relationships between ubiquitin multimers.

**Advanced Features:**
- **Graph Theory Integration:** Implements subgraph containment analysis
- **Configurable Parameters:** Allows selection of lysine types and multimer sizes
- **Matrix Visualization:** Displays containment relationships in tabular format
- **Performance Optimization:** Includes timing information and progress indicators
- **Streaming Support:** Handles long-running analysis tasks with real-time updates

### 12. `EdgeTreeViewer.jsx` - Nomenclature Visualization
**Role:** Specialized component for displaying nomenclature tree structures.

**Specialized Functions:**
- **Tree Construction:** Converts formatted edge data into hierarchical tree structures
- **Multiple Nomenclatures:** Displays all 5 nomenclature systems simultaneously
- **Subscript Formatting:** Handles complex chemical notation with subscripts
- **Data Integration:** Processes nomenclature data from backend API responses

### 13. `Panel.jsx` - UI Container Component
**Role:** Consistent styling wrapper for all major UI sections.

**Design System:**
- **Consistent Theming:** Dark theme with consistent colors (#181818 background, #444444 borders)
- **Flexible Sizing:** Supports custom dimensions while maintaining design consistency
- **Prop Forwarding:** Passes through additional props for extended functionality

---

## Design System Architecture

### Color Palette & Theming
```css
/* Dark Theme (Primary) */
:root {
  --bg-primary: #242424;
  --bg-secondary: #181818;
  --border-color: #444444;
  --text-primary: rgba(255, 255, 255, 0.87);
  --accent-blue: #646cff;
  --accent-blue-hover: #535bf2;
}

/* Light Theme (Override) */
@media (prefers-color-scheme: light) {
  :root {
    --bg-primary: #ffffff;
    --text-primary: #213547;
    --accent-blue-hover: #747bff;
  }
}
```

### Typography System
- **Primary Font Stack:** system-ui, Avenir, Helvetica, Arial, sans-serif
- **Heading Hierarchy:** h1 at 3.2em with 1.1 line-height
- **Body Text:** 1.5 line-height for optimal readability
- **Interactive Elements:** Inherit font family with 500 weight

### Layout Principles
- **Centered Container:** Max-width 1280px with auto margins
- **Consistent Spacing:** 2rem base padding, 0.6em-1.2em for buttons
- **Responsive Breakpoints:** Automatic adaptation based on content
- **Flexbox Centered:** Body uses flex with place-items: center

### Animation Standards
- **Transition Duration:** 300ms for hover effects, 0.25s for border changes
- **Performance Optimization:** will-change property for animated elements
- **Accessibility Compliance:** Respects prefers-reduced-motion settings
- **Smooth Interactions:** Hardware-accelerated transforms and filters

---

## Data Flow Architecture

### Application Hierarchy
```
main.jsx (Entry Point)
├── index.css (Global Styles)
├── App.jsx (Application Shell)
│   ├── App.css (Component Styles)
│   └── ModuleDashboard.jsx (Central State)
│       ├── Page Selection State
│       ├── API Response Data
│       ├── Interactive Scaffold State
│       └── Child Components (Specialized State)
```

### State Management Hierarchy
```
ModuleDashboard (Central State)
├── Page Selection State
├── API Response Data
├── Interactive Scaffold State
└── Selected Panels for Synthesis

Child Components (Specialized State)
├── ClickableScaffoldPanel (Canvas interactions)
├── ScaffoldPanel (Rendering state)
├── ReactionPathStatisticsPage (Chart data)
└── SubgraphAnalysisPage (Analysis results)
```

### CSS Cascade Strategy
1. **Global Reset:** index.css establishes base styles and CSS custom properties
2. **Component Overrides:** App.css provides component-specific styling
3. **Inline Styles:** Strategic use for dynamic or component-specific adjustments
4. **Component CSS:** Individual components handle their internal styling

---

## Technical Implementation Details

### Performance Optimizations
- **CSS Custom Properties:** Efficient theme switching and style management
- **Hardware Acceleration:** Transform and filter animations use GPU
- **Font Loading:** System font stack eliminates web font loading delays
- **Responsive Images:** will-change property for smooth logo animations

### Cross-Browser Compatibility
- **Font Smoothing:** -webkit-font-smoothing and -moz-osx-font-smoothing
- **Modern CSS Features:** CSS custom properties with fallbacks
- **Flexbox Layout:** Widely supported flexible layout system
- **Media Query Support:** Color scheme and motion preference detection

### Accessibility Features
- **Semantic HTML:** Proper heading hierarchy and semantic elements
- **Focus Management:** Visible focus indicators for keyboard navigation
- **Motion Preferences:** Respects user's animation preferences
- **Color Contrast:** High contrast ratios in both light and dark themes
- **Screen Reader Support:** Semantic structure and proper ARIA usage

---

## Extensibility and Maintenance

### Adding New Global Styles
1. Update CSS custom properties in :root selector
2. Add theme overrides for light mode if needed
3. Consider component-specific styles in App.css
4. Test accessibility and cross-browser compatibility

### Component Integration Guidelines
- **Follow Design System:** Use established color palette and spacing
- **Maintain Accessibility:** Include proper focus management and semantic structure
- **Performance Awareness:** Use hardware acceleration for animations
- **Theme Compatibility:** Support both light and dark modes

### Build and Development
- **Vite Configuration:** Optimized for fast development and production builds
- **Hot Module Replacement:** Instant updates during development
- **CSS Processing:** Automatic vendor prefixing and optimization
- **Asset Optimization:** Efficient bundling and code splitting

---

**This architecture provides a comprehensive, scientifically-focused interface that scales from simple visualization to complex analysis workflows, with robust styling systems, accessibility compliance, and modern development practices.**

---

## Core Architecture

### 1. `ModuleDashboard.jsx` - Main Dashboard Controller
**Role:** Central application controller and page router.

**Key Features:**
- **Page Management:** Manages 6 pages through PAGE_CONFIG: draw, tetramers, pentamers, reactionStats, subgraph, nomenclature
- **State Management:** Coordinates state across all child components including:
  - Selected panels for synthesis planning
  - API response data (figures, reaction sequences, nomenclature data)
  - Interactive scaffold state (nodes, arrows, mapping)
- **API Integration:** Handles all backend communication for synthesis planning and nomenclature analysis
- **Dynamic Rendering:** Conditionally renders different page layouts based on current page selection

**Key Imports & Dependencies:**
- `Panel` - Container styling
- `JsonToScaffold` - JSON-to-visual conversion
- `ClickableScaffoldPanel` - Interactive drawing interface
- `ReactionSequencesPaneled` - Synthesis pathway visualization
- `SubgraphAnalysisPage` - Isomorphism analysis
- `ReactionPathStatisticsPage` - Metrics visualization
- `EdgeTreeViewer` - Nomenclature tree display

---

## Page-Specific Components

### 2. `ClickableScaffoldPanel.jsx` - Interactive Drawing Interface
**Role:** Primary interactive component for the "Explore Reaction Pathways" page.

**Core Functionality:**
- **Canvas-Based Drawing:** Renders interactive node-edge graphs on HTML5 canvas
- **Click Logic:** Implements complex click-to-arrow logic for pathway construction
- **JSON Simulation:** Includes `simulateClicksFromJson()` function that converts JSON ubiquitin structures to visual representations
- **State Synchronization:** Maintains nodes, arrows, clickedNodes, and smallNodes state
- **API Integration:** Submits constructed pathways to backend for analysis

**Technical Details:**
- Uses refs for performance optimization in canvas operations
- Supports both interactive and frozen modes
- Handles complex graph traversal algorithms for pathway visualization
- Exports `simulateClicksFromJson` for use by other components

### 3. `JsonToScaffold.jsx` - Static Scaffold Renderer
**Role:** Converts JSON ubiquitin data into static visual representations.

**Key Features:**
- **JSON Processing:** Implements its own `simulateClicksFromJson()` function (similar to ClickableScaffoldPanel)
- **Preorder Traversal:** Follows biological conventions for ubiquitin chain numbering
- **Visual Logic:** K63 linkages move rightward/up, K48 linkages move leftward/up
- **Protection Groups:** Renders SMAC/ABOC protecting groups as small colored nodes
- **Static Display:** Uses ScaffoldPanel in frozen mode for non-interactive display

**Usage Context:**
- Used in tetramer/pentamer selection grids
- Renders individual multimer thumbnails
- Provides consistent visual representation across the application

### 4. `ScaffoldPanel.jsx` - Core Visualization Engine
**Role:** Low-level canvas rendering engine used by both interactive and static components.

**Technical Implementation:**
- **Canvas Management:** Handles all drawing operations (nodes, edges, arrows, small nodes)
- **State Management:** Tracks node states, arrow configurations, and small node positions
- **Event Handling:** Manages mouse interactions when not frozen
- **Visual Consistency:** Maintains consistent styling and colors across all scaffold displays
- **Modular Design:** Supports both interactive (ClickableScaffoldPanel) and static (JsonToScaffold) use cases

### 5. `ReactionSequencesPaneled.jsx` - Synthesis Pathway Display
**Role:** Visualizes step-by-step reaction sequences for synthesis planning.

**Features:**
- **Multi-Step Visualization:** Shows progression from dimer → trimer → tetramer → pentamer
- **Reaction Well Display:** Optional display of specific reaction conditions
- **JSON Processing:** Handles complex nested JSON structures from backend
- **Dynamic Layout:** Adapts to both tetramer and pentamer synthesis pathways
- **Integration:** Uses JsonToScaffold for rendering individual synthesis steps

### 6. `ReactionPathStatisticsPage.jsx` - Metrics Dashboard
**Role:** Statistical analysis and visualization for reaction pathways.

**Capabilities:**
- **Chart Integration:** Uses Chart.js for bar charts and scatter plots
- **API Communication:** Fetches statistical data from backend endpoints
- **Interactive Controls:** Allows filtering by multimer size (4 or 5) and pathway type (ABOC/SMAC)
- **Data Processing:** Handles complex statistical datasets and presents them visually

### 7. `SubgraphAnalysisPage.jsx` - Isomorphism Analysis
**Role:** Advanced analysis of structural relationships between ubiquitin multimers.

**Advanced Features:**
- **Graph Theory Integration:** Implements subgraph containment analysis
- **Configurable Parameters:** Allows selection of lysine types and multimer sizes
- **Matrix Visualization:** Displays containment relationships in tabular format
- **Performance Optimization:** Includes timing information and progress indicators
- **Streaming Support:** Handles long-running analysis tasks with real-time updates

### 8. `EdgeTreeViewer.jsx` - Nomenclature Visualization
**Role:** Specialized component for displaying nomenclature tree structures.

**Specialized Functions:**
- **Tree Construction:** Converts formatted edge data into hierarchical tree structures
- **Multiple Nomenclatures:** Displays all 5 nomenclature systems simultaneously
- **Subscript Formatting:** Handles complex chemical notation with subscripts
- **Data Integration:** Processes nomenclature data from backend API responses

### 9. `Panel.jsx` - UI Container Component
**Role:** Consistent styling wrapper for all major UI sections.

**Design System:**
- **Consistent Theming:** Dark theme with consistent colors (#181818 background, #444444 borders)
- **Flexible Sizing:** Supports custom dimensions while maintaining design consistency
- **Prop Forwarding:** Passes through additional props for extended functionality

---

## Data Flow Architecture

### State Management Hierarchy
```
ModuleDashboard (Central State)
├── Page Selection State
├── API Response Data
├── Interactive Scaffold State
└── Selected Panels for Synthesis

Child Components (Specialized State)
├── ClickableScaffoldPanel (Canvas interactions)
├── ScaffoldPanel (Rendering state)
├── ReactionPathStatisticsPage (Chart data)
└── SubgraphAnalysisPage (Analysis results)
```

### API Integration Points
1. **Synthesis Planning:** Submit selected tetramers/pentamers → receive synthesis protocols
2. **Pathway Exploration:** Submit scaffold configurations → receive reaction sequences
3. **Nomenclature Analysis:** Submit UbID/multimer requests → receive nomenclature data
4. **Statistical Analysis:** Request pathway metrics → receive chart data
5. **Isomorphism Analysis:** Submit analysis parameters → receive containment matrices

### Component Communication Patterns
- **Props Down:** Configuration and data flow from ModuleDashboard to specialized components
- **Callbacks Up:** User interactions and results flow back to ModuleDashboard via callback props
- **Shared Utilities:** `simulateClicksFromJson` implemented in multiple components for consistency
- **State Synchronization:** ModuleDashboard coordinates state between related components

---

## Technical Implementation Details

### JSON Processing Pipeline
1. **Input:** Backend provides nested JSON structures representing ubiquitin chains
2. **Processing:** `simulateClicksFromJson` functions convert JSON to visual coordinates
3. **Rendering:** ScaffoldPanel renders the processed data on HTML5 canvas
4. **Interaction:** ClickableScaffoldPanel handles user modifications
5. **Output:** Modified structures sent back to backend for analysis

### Visual Consistency System
- **Shared Constants:** DEFAULT_NODES, DEFAULT_EDGES defined across components
- **Color Coding:** K63 (blue), K48 (red), protection groups (colored small nodes)
- **Positioning Logic:** Biological conventions enforced in positioning algorithms
- **Responsive Design:** Components adapt to different panel sizes while maintaining proportions

### Performance Optimizations
- **Canvas Refs:** Direct DOM manipulation for smooth interactions
- **State Memoization:** React.useMemo for expensive calculations
- **Conditional Rendering:** Components only render when data is available
- **Background Processing:** Long-running analyses handled asynchronously

---

## Extensibility and Maintenance

### Adding New Pages
1. Add page configuration to PAGE_CONFIG in ModuleDashboard
2. Create new component following established patterns
3. Integrate API endpoints if needed
4. Add navigation option to page selector

### Component Modification Guidelines
- **Scaffold Components:** Maintain consistency in visual representation
- **API Components:** Follow established error handling patterns
- **UI Components:** Respect the dark theme design system
- **State Management:** Use callback patterns for parent-child communication

### Testing Considerations
- **Canvas Testing:** Visual components require specialized testing approaches
- **API Mocking:** Backend dependencies should be mockable for unit tests
- **State Testing:** Complex state interactions benefit from integration tests
- **Performance Testing:** Canvas operations should be performance tested

---

**This architecture provides a robust, scientifically-focused interface that scales from simple visualization to complex analysis workflows, with clear separation of concerns and consistent user experience patterns.**
