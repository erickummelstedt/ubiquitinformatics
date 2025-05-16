import { render, screen } from '@testing-library/react';
import SequenceFilter from '../../src/components/SequenceFilter';

test('renders filter input', () => {
  render(<SequenceFilter />);
  expect(screen.getByPlaceholderText(/Filter sequences/i)).toBeInTheDocument();
});
