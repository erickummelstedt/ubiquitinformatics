import { render, screen } from '@testing-library/react';
import SequenceList from '../../src/components/SequenceList';

test('renders sequence list', () => {
  const sequences = [{ name: 'Test Sequence' }];
  render(<SequenceList sequences={sequences} onSelect={() => {}} />);
  expect(screen.getByText(/Test Sequence/i)).toBeInTheDocument();
});
