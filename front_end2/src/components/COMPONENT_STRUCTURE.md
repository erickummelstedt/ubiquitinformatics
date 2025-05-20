# Component Structure and Connections

## Overview
This front-end React/Vite dashboard is organized into modular components for clarity, reusability, and maintainability. Below is a description of each main component and how they connect to form the interactive dashboard for ubiquitin multimer synthesis simulation.

---

## 1. `SequenceFilter.jsx`
**Role:** Main layout and composition file.

- **Imports:**
  - `Panel` (panel styling/layout)
  - `GameScaffoldPanel` (interactive scaffold/canvas)
  - `Visualizer` (3Dmol.js molecular viewer)
- **Layout:**
  - Uses a flex row to display two panels side by side, with a gap for visual separation.
  - **Left Panel:** Contains the interactive scaffold (GameScaffoldPanel).
  - **Right Panel:** Contains the molecular visualizer (Visualizer).
- **Props:**
  - Passes the .mol2 file URL as a prop to `Visualizer`.

---

## 2. `Panel.jsx`
**Role:** Generic styled container for content panels.

- **Usage:** Wraps both the left and right panels in `SequenceFilter`.
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
  - Used in the left panel of `SequenceFilter`.

---

## 4. `Visualizer.jsx`
**Role:** 3D molecular structure viewer using 3Dmol.js.

- **Features:**
  - Loads and displays a .mol2 file in a contained 3D viewer.
  - Strictly contained within its parent panel.
- **Props:**
  - `mol2Url`: Path/URL to the .mol2 file to visualize.
- **Logic:**
  - Dynamically loads 3Dmol.js if needed.
  - Fetches and renders the molecule.
- **Usage:**
  - Used in the right panel of `SequenceFilter`.

---

## How It All Connects
- `SequenceFilter.jsx` is the main entry point for the dashboard UI.
- It composes two `Panel` components side by side:
  - The left panel contains `GameScaffoldPanel` for interactive simulation.
  - The right panel contains `Visualizer` for molecular structure viewing.
- All components are modular and can be edited or extended independently.
- Data and configuration (e.g., .mol2 file path, scaffold structure) are passed as props for flexibility.

---

## Extensibility
- New panels or features can be added by creating new components and including them in `SequenceFilter`.
- Logic can be further abstracted into custom hooks for testability.
- Styling can be moved to CSS modules or styled-components for theming.

---

**This structure ensures a robust, maintainable, and extensible codebase for interactive scientific dashboards.**
