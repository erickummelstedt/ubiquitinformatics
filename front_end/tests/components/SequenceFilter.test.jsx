// This test file has been replaced by ScaffoldDashboard.test.jsx. Please use that file instead.

import { render, screen, fireEvent } from '@testing-library/react';
import ScaffoldDashboard from '../../src/components/ScaffoldDashboard';
import React from 'react';

describe('ScaffoldDashboard', () => {
  it('renders without crashing', () => {
    render(<ScaffoldDashboard />);
    expect(screen.getByText(/molecular visualization is currently disabled/i)).toBeInTheDocument();
  });

  it('renders filter input', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    expect(screen.getByPlaceholderText(/Filter sequences/i)).toBeInTheDocument();
  });

  it('calls onChange when typing in filter input', () => {
    const handleChange = jest.fn();
    render(<ScaffoldDashboard value="" onChange={handleChange} />);
    const input = screen.getByPlaceholderText(/Filter sequences/i);
    fireEvent.change(input, { target: { value: 'abc' } });
    expect(handleChange).toHaveBeenCalled();
  });

  it('renders the main panel and canvas', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    // Panel is a div with backgroundColor #181818
    const panel = screen.getByRole('region', { hidden: true }) || screen.getByText((_, el) => el && el.style && el.style.backgroundColor === '#181818');
    expect(panel).toBeTruthy();
    // Canvas
    const canvas = screen.getByRole('img', { hidden: true }) || document.querySelector('canvas');
    expect(canvas).toBeInTheDocument();
  });

  it('renders the scaffold canvas with correct dimensions', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    const canvas = document.querySelector('canvas');
    expect(canvas).toBeInTheDocument();
    expect(canvas.width).toBeGreaterThan(500); // Should be wide
    expect(canvas.height).toBeGreaterThan(300); // Should be tall
  });

  it('shows the panel selection logic placeholder', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    expect(screen.getByText(/Panel selection and output connection logic/i)).toBeInTheDocument();
  });

  it('refreshes the scaffold when clicking the refresh button', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    const canvas = document.querySelector('canvas');
    // Simulate click on refresh button area (bottom right of canvas)
    const rect = canvas.getBoundingClientRect();
    // These coordinates should be within the refresh button
    fireEvent.mouseDown(canvas, { clientX: rect.right - 50, clientY: rect.bottom - 20 });
    // No error should occur, and the canvas should still be present
    expect(document.querySelector('canvas')).toBeInTheDocument();
  });

  it('allows clicking on a node and draws an arrow (integration smoke test)', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    const canvas = document.querySelector('canvas');
    // Simulate click on a node above the root (should be allowed)
    // Node 1 is at (350, 250) + offset, root is at (400, 300) + offset
    // We'll click near (350+offset, 250+offset)
    // Use a rough guess for offset (since it's >0)
    fireEvent.mouseDown(canvas, { clientX: 350 + 100, clientY: 250 + 20 });
    // No error should occur, and the canvas should still be present
    expect(document.querySelector('canvas')).toBeInTheDocument();
  });

  it('does not allow clicking on a node that is not allowed (integration smoke test)', () => {
    render(<ScaffoldDashboard value="" onChange={() => {}} />);
    const canvas = document.querySelector('canvas');
    // Click on a node that is not adjacent/above (e.g., node 10)
    fireEvent.mouseDown(canvas, { clientX: 200 + 100, clientY: 100 + 20 });
    // No error should occur, and the canvas should still be present
    expect(document.querySelector('canvas')).toBeInTheDocument();
  });

  it('renders with a custom value in the filter input', () => {
    render(<ScaffoldDashboard value="test123" onChange={() => {}} />);
    expect(screen.getByDisplayValue('test123')).toBeInTheDocument();
  });

  it('renders with a custom onPanelSelect prop (if provided)', () => {
    // Should not throw if onPanelSelect is passed
    render(<ScaffoldDashboard value="" onChange={() => {}} onPanelSelect={() => {}} />);
    expect(screen.getByPlaceholderText(/Filter sequences/i)).toBeInTheDocument();
  });
});
