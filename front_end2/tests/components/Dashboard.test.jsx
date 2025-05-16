import { render, screen } from '@testing-library/react';
import Dashboard from '../../src/components/Dashboard';

test('renders dashboard title', () => {
  render(<Dashboard />);
  expect(screen.getByText(/Ubiquitin Multimer Synthesis Dashboard/i)).toBeInTheDocument();
});
