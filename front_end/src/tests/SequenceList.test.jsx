import React from 'react';
import { render, screen } from '@testing-library/react';
import SequenceList from './SequenceList.jsx';

describe('SequenceList', () => {
  const sequences = [
    { id: 1, name: 'Sequence 1' },
    { id: 2, name: 'Sequence 2' },
  ];

  test('renders SequenceList component', () => {
    render(<SequenceList sequences={sequences} />);

    const sequenceElements = screen.getAllByRole('listitem');
    expect(sequenceElements).toHaveLength(sequences.length);
  });

  test('displays the correct sequence names', () => {
    render(<SequenceList sequences={sequences} />);

    sequences.forEach((sequence) => {
      const sequenceName = screen.getByText(sequence.name);
      expect(sequenceName).toBeInTheDocument();
    });
  });

  test('renders empty state when no sequences are present', () => {
    render(<SequenceList sequences={[]} />);

    const emptyStateMessage = screen.getByText(/no sequences found/i);
    expect(emptyStateMessage).toBeInTheDocument();
  });
});