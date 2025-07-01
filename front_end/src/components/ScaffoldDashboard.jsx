import React, { useState } from 'react';
import Panel from './Panel';
import GameScaffoldPanel from './GameScaffoldPanel';
import FrozenGameScaffoldPanel from './FrozenGameScaffoldPanel';

const SMALL_PANEL_WIDTH = 140;
const SMALL_PANEL_HEIGHT = 90;

const PAGE_CONFIG = {
  draw: { count: 1, label: 'Explore Ubiquitin Pathways', panelWidth: 570, panelHeight: 370 },
  tetramers: { count: 14, label: 'Tetramers', panelWidth: SMALL_PANEL_WIDTH, panelHeight: SMALL_PANEL_HEIGHT },
  pentamers: { count: 42, label: 'Pentamers', panelWidth: SMALL_PANEL_WIDTH, panelHeight: SMALL_PANEL_HEIGHT },
};

const ScaffoldDashboard = () => {
  const [page, setPage] = useState('draw');
  const [selectedPanels, setSelectedPanels] = useState([]); // Track selected panel indices
  const { count, panelWidth, panelHeight } = PAGE_CONFIG[page];

  // Reset selection when page changes
  React.useEffect(() => {
    setSelectedPanels([]);
  }, [page]);

  // Count selections for each label
  const selectionCounts = {};
  selectedPanels.forEach(idx => {
    const label = page === 'tetramers' ? `Ub4_${idx + 1}` : page === 'pentamers' ? `Ub5_${idx + 1}` : null;
    if (label) selectionCounts[label] = (selectionCounts[label] || 0) + 1;
  });

  const handlePanelClick = i => {
    if (page === 'draw') return;
    if (selectedPanels.length >= 16) return; // Limit to 16 selections
    setSelectedPanels(prev => [...prev, i]);
  };

  // Submit handler to call FastAPI backend
  const handleSubmit = async () => {
    // Prepare the selected label list (e.g., ["Ub4_1", "Ub4_2", ...])
    const selectedLabels = selectedPanels.map(idx =>
      page === 'tetramers' ? `Ub4_${idx + 1}` : page === 'pentamers' ? `Ub5_${idx + 1}` : null
    ).filter(Boolean);
    try {
      const response = await fetch('/api/submit-selection', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ labels: selectedLabels, page }),
      });
      if (!response.ok) throw new Error('Submission failed');
      const result = await response.json();
      alert('Submission successful!\n' + JSON.stringify(result));
    } catch (err) {
      alert('Submission failed: ' + err.message);
    }
  };

  return (
    <div style={{ width: '100vw', minHeight: '100vh', boxSizing: 'border-box', overflow: 'auto', padding: '16px 0' }}>
      <div style={{ display: 'flex', justifyContent: 'center', marginBottom: 24 }}>
        <label htmlFor="dashboard-page-select" style={{ fontWeight: 600, fontSize: 16, marginRight: 8 }}>Page:</label>
        <select
          id="dashboard-page-select"
          value={page}
          onChange={e => setPage(e.target.value)}
          style={{ fontSize: 16, padding: '4px 12px', borderRadius: 6 }}
        >
          <option value="draw">Draw</option>
          <option value="tetramers">Tetramers</option>
          <option value="pentamers">Pentamers</option>
        </select>
      </div>
      <div style={{
        display: 'flex',
        flexWrap: 'wrap',
        alignItems: 'flex-start',
        justifyContent: 'center',
        gap: page === 'draw' ? '0' : '8px',
        rowGap: page === 'draw' ? '0' : '8px',
        columnGap: page === 'draw' ? '0' : '8px',
      }}>
        {[...Array(count)].map((_, i) => {
          const isSelected = selectedPanels.includes(i);
          const label = page === 'tetramers' ? `Ub4_${i + 1}` : page === 'pentamers' ? `Ub5_${i + 1}` : '';
          return (
            <Panel
              key={i}
              style={{
                width: panelWidth,
                height: panelHeight + (page === 'draw' ? 0 : 20),
                minWidth: panelWidth,
                minHeight: panelHeight + (page === 'draw' ? 0 : 20),
                padding: 0,
                position: 'relative',
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                margin: page === 'draw' ? '0 auto' : undefined,
                border: isSelected ? '3px solid #1976d2' : '2px solid transparent',
                boxShadow: isSelected ? '0 0 12px #1976d2' : undefined,
                cursor: page !== 'draw' ? 'pointer' : 'default',
                transition: 'border 0.2s, box-shadow 0.2s',
              }}
              onClick={() => handlePanelClick(i)}
            >
              <div style={{ width: panelWidth, height: panelHeight }}>
                {page === 'draw' ? (
                  <GameScaffoldPanel panelWidth={panelWidth} panelHeight={panelHeight} />
                ) : (
                  <FrozenGameScaffoldPanel panelWidth={panelWidth} panelHeight={panelHeight} />
                )}
              </div>
              {page !== 'draw' && (
                <div
                  style={{
                    width: '100%',
                    textAlign: 'center',
                    fontSize: 14,
                    color: 'white',
                    marginTop: -5,
                    fontWeight: 600,
                    textShadow: '0 2px 8px #222',
                    letterSpacing: 1,
                    userSelect: 'none',
                    position: 'relative',
                    zIndex: 2,
                  }}
                >
                  {label}
                  {selectionCounts[label] ? (
                    <span style={{ marginLeft: 8, color: '#ffd600', fontWeight: 700, fontSize: 13 }}>
                      ×{selectionCounts[label]}
                    </span>
                  ) : null}
                </div>
              )}
            </Panel>
          );
        })}
      </div>
      {/* Selection tracker for tetramers and pentamers */}
      {page !== 'draw' && (
        <div style={{
          width: '100%',
          maxWidth: 900,
          margin: '32px auto 0 auto',
          background: '#181818',
          borderRadius: 10,
          padding: '16px 24px',
          color: 'white',
          fontSize: 16,
          fontWeight: 500,
          boxShadow: '0 2px 12px #0006',
          minHeight: 40,
          display: 'flex',
          alignItems: 'center',
          gap: 16,
        }}>
          <span style={{ fontWeight: 700, marginRight: 12 }}>Selected:</span>
          {Object.keys(selectionCounts).length === 0 ? (
            <span style={{ color: '#aaa' }}>None</span>
          ) : (
            Object.entries(selectionCounts).map(([label, count], idx) => (
              <span key={label} style={{
                display: 'inline-block',
                background: '#1976d2',
                color: '#fff',
                borderRadius: 6,
                padding: '2px 10px',
                marginRight: 8,
                fontWeight: 700,
                fontSize: 15,
                boxShadow: '0 1px 4px #0004',
              }}>
                {label}{count > 1 ? ` ×${count}` : ''}
              </span>
            ))
          )}
          <button
            onClick={() => setSelectedPanels([])}
            style={{
              marginLeft: 'auto',
              background: '#fff',
              color: '#1976d2',
              border: 'none',
              borderRadius: 6,
              padding: '6px 18px',
              fontWeight: 700,
              fontSize: 15,
              cursor: 'pointer',
              boxShadow: '0 1px 4px #0002',
              transition: 'background 0.2s, color 0.2s',
            }}
            title="Clear all selections"
          >
            Refresh
          </button>
          <button
            onClick={handleSubmit}
            style={{
              marginLeft: 12,
              background: '#1976d2',
              color: '#fff',
              border: 'none',
              borderRadius: 6,
              padding: '6px 18px',
              fontWeight: 700,
              fontSize: 15,
              cursor: 'pointer',
              boxShadow: '0 1px 4px #0002',
              transition: 'background 0.2s, color 0.2s',
            }}
            disabled={selectedPanels.length === 0}
            title="Submit selection to backend"
          >
            Submit
          </button>
        </div>
      )}
    </div>
  );
};

export default ScaffoldDashboard;
