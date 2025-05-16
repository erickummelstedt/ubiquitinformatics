import { render, screen } from '@testing-library/react';
import Visualizer from '../../src/components/Visualizer';

test('renders visualizer placeholder', () => {
  render(<Visualizer data={null} />);
  expect(screen.getByText(/Visualizer placeholder/i)).toBeInTheDocument();
});
