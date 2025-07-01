# Component Structure and Connections

## Overview
This front-end React/Vite dashboard is organized into modular components for clarity, reusability, and maintainability. Below is a description of each main component and how they connect to form the interactive dashboard for ubiquitin multimer synthesis simulation.

---

## 1. `ScaffoldDashboard.jsx`
**Role:** Main layout and composition file.

- **Imports:**
  - `Panel` (panel styling/layout)
  - `GameScaffoldPanel` (interactive scaffold/canvas)
- **Layout:**
  - Uses a flex row to display a single panel for the interactive scaffold.
  - **Panel:** Contains the interactive scaffold (`GameScaffoldPanel`).
- **Props:**
  - None specific (can be extended as needed).

---

## 2. `Panel.jsx`
**Role:** Generic styled container for content panels.

- **Usage:** Wraps the main panel in `ScaffoldDashboard`.
- **Props:**
  - `children`: Content to display inside the panel.
  - `style`: Optional style overrides.
- **Styling:**
  - Consistent background, border, padding, and rounded corners.
  - Fixed width and height for layout consistency.

---

## 3. `GameScaffoldPanel.jsx`
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
  - Used in the main panel of `ScaffoldDashboard`.

---

## How It All Connects
- `ScaffoldDashboard.jsx` is the main entry point for the dashboard UI.
- It composes a single `Panel` component containing `GameScaffoldPanel` for interactive simulation.
- All components are modular and can be edited or extended independently.
- Data and configuration (e.g., scaffold structure) are passed as props for flexibility.

---

## Extensibility
- New panels or features can be added by creating new components and including them in `ScaffoldDashboard`.
- Logic can be further abstracted into custom hooks for testability.
- Styling can be moved to CSS modules or styled-components for theming.

---

**This structure ensures a robust, maintainable, and extensible codebase for interactive scientific dashboards.**
