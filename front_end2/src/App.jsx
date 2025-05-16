import React, { useEffect, useState } from 'react';
import SequenceFilter from './components/SequenceFilter';
import SequenceList from './components/SequenceList';
import Visualizer from './components/Visualizer';
import { fetchInitialReactions } from './services/api';
import './index.css';

function App() {
  const [sequences, setSequences] = useState([]);
  const [selected, setSelected] = useState(null);
  const [filter, setFilter] = useState('');

  useEffect(() => {
    fetchInitialReactions().then(data => {
      setSequences(data);
      setSelected(data[0] || null);
    });
  }, []);

  // Simple filter for demonstration (can be improved)
  const filteredSequences = sequences.filter(seq =>
    Object.values(seq).some(val => val.toLowerCase().includes(filter.toLowerCase()))
  );

  return (
    <div className="App">
      <h1>Ubiquitin Multimer Synthesis Dashboard</h1>
      <SequenceFilter value={filter} onChange={e => setFilter(e.target.value)} />
      <SequenceList sequences={filteredSequences} onSelect={setSelected} />
      <Visualizer data={selected} />
    </div>
  );
}

export default App;
