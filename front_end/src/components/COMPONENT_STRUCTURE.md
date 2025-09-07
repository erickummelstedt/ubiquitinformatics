# Component Structure and Connections

## Overview
This front-end React/Vite dashboard is organized into modular components for clarity, reusability, and maintainability. Below is a description of each main component and how they connect to form the interactive dashboard for ubiquitin multimer synthesis simulation.

---

## 1. `ModuleDashboard.jsx`
**Role:** Main layout and composition file.

- **Imports:**
  - `Panel` (panel styling/layout)
  - `ClickableScaffoldPanel` (interactive scaffold/canvas)
  - `JsonToScaffold` (JSON-based scaffold rendering)
  - `Sequences` (sequence visualization)
  - `SubgraphAnalysisPage` (subgraph analysis UI)
  - `ReactionPathStatisticsPage` (reaction path metrics visualization)
- **Layout:**
  - Dynamically renders panels based on the selected page.
  - **Panels:** Includes interactive scaffold, reaction metrics, and subgraph analysis.
- **Props:**
  - None specific (can be extended as needed).

---

## 2. `Panel.jsx`
**Role:** Generic styled container for content panels.

- **Usage:** Wraps the main panel in `ModuleDashboard`.
- **Props:**
  - `children`: Content to display inside the panel.
  - `style`: Optional style overrides.
- **Styling:**
  - Consistent background, border, padding, and rounded corners.
  - Fixed width and height for layout consistency.

---

## 3. `ClickableScaffoldPanel.jsx`
**Role:** Interactive, game-like scaffold visualization using a canvas.

- **Features:**
  - Renders a tree/graph of nodes and edges representing a synthesis pathway.
  - Nodes are clickable; clicking adds arrows (red/blue) according to logic.
  - Includes a refresh button to reset the scaffold.
- **Props:**
  - `initialNodes` (optional): Array of node objects (default provided).
  - `edges` (optional): Array of edge pairs (default provided).
- **State:**
  - Tracks node states, arrows, clicked nodes, and refresh button.
- **Logic:**
  - Handles drawing, click detection, and interactive updates.
- **Usage:**
  - Used in the main panel of `ModuleDashboard`.

---

## 4. `JsonToScaffold.jsx`
**Role:** Converts JSON data into scaffold visualizations.

- **Features:**
  - Simulates clicks and renders scaffold structures based on JSON input.
- **Props:**
  - `jsonData`: JSON object representing scaffold data.
- **Usage:**
  - Used in `ModuleDashboard` and other components for scaffold rendering.

---

## 5. `ReactionPathStatisticsPage.jsx`
**Role:** Displays reaction path metrics using charts.

- **Features:**
  - Bar charts and scatter plots for reaction metrics.
  - Fetches data from the backend API.
- **Props:**
  - None specific.
- **Usage:**
  - Rendered in `ModuleDashboard` for reaction path analysis.

---

## 6. `ScaffoldPanel.jsx`
**Role:** Static scaffold visualization with optional interactivity.

- **Features:**
  - Renders nodes, edges, and arrows.
  - Supports frozen state for static display.
- **Props:**
  - `initialNodes`, `edges`, `arrows`: Define scaffold structure.
  - `frozen`: Disables interactivity.
- **Usage:**
  - Used in `JsonToScaffold` and other scaffold-related components.

---

## 7. `Sequences.jsx`
**Role:** Renders detailed reaction pathways and scaffolds.

- **Features:**
  - Displays reaction wells and multimer synthesis steps.
  - Uses `JsonToScaffold` for scaffold rendering.
- **Props:**
  - `reactionSequence`: Array of reaction data.
  - `showReactionWell`: Toggles reaction well display.
- **Usage:**
  - Used in `ModuleDashboard` for detailed pathway visualization.

---

## 8. `SubgraphAnalysisPage.jsx`
**Role:** Provides UI for subgraph containment analysis.

- **Features:**
  - Configurable lysine IDs and multimer sizes.
  - Displays analysis results in a table.
- **Props:**
  - None specific.
- **Usage:**
  - Rendered in `ModuleDashboard` for subgraph analysis.

---

## How It All Connects
- `ModuleDashboard.jsx` is the main entry point for the dashboard UI.
- It dynamically renders panels based on the selected page.
- Components like `JsonToScaffold`, `ReactionPathStatisticsPage`, and `SubgraphAnalysisPage` provide specialized functionality.
- All components are modular and can be edited or extended independently.
- Data and configuration (e.g., scaffold structure) are passed as props for flexibility.

---

## Extensibility
- New panels or features can be added by creating new components and including them in `ModuleDashboard`.
- Logic can be further abstracted into custom hooks for testability.
- Styling can be moved to CSS modules or styled-components for theming.

---

**This structure ensures a robust, maintainable, and extensible codebase for interactive scientific dashboards.**
