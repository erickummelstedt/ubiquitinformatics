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
      // If the data is an array of objects, convert each object to an array of {key, value} pairs
      const processed = data.map(obj =>
        Object.entries(obj).map(([key, value]) => ({ key, value }))
      );
      setSequences(processed);
      setSelected(processed[0] || null);
    });
  }, []);

  // Simple filter for demonstration (can be improved)
  const filteredSequences = sequences.filter(seqArr =>
    seqArr.some(({ key, value }) =>
      key.toLowerCase().includes(filter.toLowerCase()) ||
      String(value).toLowerCase().includes(filter.toLowerCase())
    )
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
