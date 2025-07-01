import React from 'react';
import { render, screen } from '@testing-library/react';
import ScaffoldDashboard from '../../src/components/ScaffoldDashboard';

describe('ScaffoldDashboard', () => {
  it('renders without crashing', () => {
    render(<ScaffoldDashboard />);
    expect(screen.getByText(/molecular visualization is currently disabled/i)).toBeInTheDocument();
  });
});
