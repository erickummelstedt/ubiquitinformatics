import React from 'react';
import ScaffoldDashboard from './components/ScaffoldDashboard';
import './index.css';

function App() {
  return (
    <div className="App">
      <h1 style={{ textAlign: 'center', margin: '32px 0 24px 0', color: 'black', fontWeight: 'bold', textShadow: '0 2px 8px #fff' }}>
        Ubiquitin Multimer Synthesis Dashboard
      </h1>
      <ScaffoldDashboard />
    </div>
  );
}

export default App;
