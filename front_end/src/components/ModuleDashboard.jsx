import React, { useState, useEffect } from 'react';
import Panel from './Panel';
import JsonToScaffold from './JsonToScaffold';
import ClickableScaffoldPanel, { simulateClicksFromJson } from './ClickableScaffoldPanel';
import ReactionSequencesPaneled from './ReactionSequencesPaneled';
import SubgraphAnalysisPage from './SubgraphAnalysisPage';
import ReactionPathStatisticsPage from './ReactionPathStatisticsPage';
import EdgeTreeViewer from './EdgeTreeViewer';
import multimerDataTetramers from '../data/multimer_id_to_json4.json';
import multimerDataPentamers from '../data/multimer_id_to_json5.json';

const SMALL_PANEL_WIDTH = 140;
const SMALL_PANEL_HEIGHT = 90;

const DEFAULT_NODES = [
  { x: 300, y: 300, clicks: 0 },
  { x: 250, y: 250, clicks: 0 }, { x: 350, y: 250, clicks: 0 },
  { x: 200, y: 200, clicks: 0 }, { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 },
  { x: 150, y: 150, clicks: 0 }, { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 },
  { x: 100, y: 100, clicks: 0 }, { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }
];

const TEST_NODES = [
  { x: 300, y: 300, clicks: 1 },
  { x: 250, y: 250, clicks: 0 }, { x: 350, y: 250, clicks: 0 },
  { x: 200, y: 200, clicks: 0 }, { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 },
  { x: 150, y: 150, clicks: 1 }, { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 },
  { x: 100, y: 100, clicks: 0 }, { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }
];

const DEFAULT_EDGES = [
  [0, 1], [0, 2],
  [1, 3], [1, 4],
  [2, 4], [2, 5],
  [3, 6], [3, 7],
  [4, 7], [4, 8],
  [5, 8], [5, 9],
  [6, 10], [6, 11],
  [7, 11], [7, 12],
  [8, 12], [8, 13],
  [9, 13], [9, 14]
];

const PAGE_CONFIG = {
  draw: { count: 1, label: 'Explore Reaction Pathways', panelWidth: 570, panelHeight: 370 },
  tetramers: { count: 14, label: 'Tetramers', panelWidth: SMALL_PANEL_WIDTH, panelHeight: SMALL_PANEL_HEIGHT },
  pentamers: { count: 42, label: 'Pentamers', panelWidth: SMALL_PANEL_WIDTH, panelHeight: SMALL_PANEL_HEIGHT },
  reactionStats: { count: 1, label: 'Reaction Path Metrics', panelWidth: 570, panelHeight: 500 },
  subgraph: { count: 1, label: 'Ubiquitin Isomorphism', panelWidth: 570, panelHeight: 500 },
};

const ModuleDashboard = () => {
  const [page, setPage] = useState('draw');
  const [selectedPanels, setSelectedPanels] = useState([]); // Track selected panel indices
  const [figures, setFigures] = useState(null); // Store backend images
  const [reactionSequence, setReactionSequence] = useState(null); // Store reaction sequences
  const [jsonOutput, setJsonOutput] = useState(null);
  const [formattedEdges, setFormattedEdges] = useState(null); // Store formatted edges from API
  const [ubxyValue, setUbxyValue] = useState(null); // Store UbX_Y value from API
  const [nomenclaturePreorderABC, setNomenclaturePreorderABC] = useState(null); // Store nomenclature value from API
  const [stritarNomenclatureWoPreorder, setStritarNomenclatureWoPreorder] = useState(null); // Store stritar nomenclature (without preorder)
  const [kummelstedtNomenclatureWPreorder, setKummelstedtNomenclatureWPreorder] = useState(null); // Store kummelstedt nomenclature (with preorder)
  const [inputNodes, setInputNodes] = useState({
    nodes: DEFAULT_NODES,
    arrows: [],
    mapping: {},
    smallNodes: []
  });
  const [scaffoldKey, setScaffoldKey] = useState(Date.now());

  const { count, panelWidth, panelHeight } = PAGE_CONFIG[page];

  // Reset selection and images when page changes
  React.useEffect(() => {
    setSelectedPanels([]);
    setFigures(null);
    setFormattedEdges(null);
    setUbxyValue(null);
    setNomenclaturePreorderABC(null);
    setStritarNomenclatureWoPreorder(null);
    setKummelstedtNomenclatureWPreorder(null);
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
    const selectedLabels = selectedPanels.map(idx =>
      page === 'tetramers' ? `Ub4_${idx + 1}` : page === 'pentamers' ? `Ub5_${idx + 1}` : null
    ).filter(Boolean);
    try {
      // Clear previous reaction sequences and scaffolds
      setReactionSequence(null);
      setFigures(null);
      setJsonOutput(null); // Clear previously rendered scaffold
      setFormattedEdges(null); // Clear formatted edges
      setUbxyValue(null); // Clear UbX_Y value
      setNomenclaturePreorderABC(null); // Clear nomenclature value
      setStritarNomenclatureWoPreorder(null); // Clear stritar nomenclature
      setKummelstedtNomenclatureWPreorder(null); // Clear kummelstedt nomenclature

      const response = await fetch('/api/submit-selection', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ labels: selectedLabels, page }),
      });
      if (!response.ok) throw new Error('Submission failed');
      const result = await response.json();
      setFigures(result.figures); // Store images
      // alert('Submission successful!\n' + JSON.stringify(result));
    } catch (err) {
      alert('Submission failed: ' + err.message);
    }
  };

  const fixQuotes = (str) => str.replace(/'/g, '"'); // Replace single quotes with double quotes

  const getMultimerData = (label) => {
    let data;
    if (page === 'tetramers') {
      data = multimerDataTetramers[label];
    } else if (page === 'pentamers') {
      data = multimerDataPentamers[label];
    }
    return data ? JSON.parse(fixQuotes(data)) : null;
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
          <option value="draw">Explore Reaction Pathways</option>
          <option value="tetramers">Tetramers</option>
          <option value="pentamers">Pentamers</option>
          <option value="reactionStats">Reaction Path Metrics</option>
          <option value="subgraph">Ubiquitin Isomorphism</option>
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
        {page === 'reactionStats' ? (
          <div style={{ 
            width: '100%', 
            padding: '0 32px', 
            boxSizing: 'border-box',
            maxWidth: '1400px'
          }}>
            <ReactionPathStatisticsPage />
          </div>
        ) : page === 'subgraph' ? (
          <div style={{ 
            width: '100%', 
            padding: '0 32px', 
            boxSizing: 'border-box',
            maxWidth: '1400px'
          }}>
            <SubgraphAnalysisPage />
          </div>
        ) : (
          [...Array(count)].map((_, i) => {
            const isSelected = selectedPanels.includes(i);
            const label = page === 'tetramers' ? `Ub4_${i + 1}` : page === 'pentamers' ? `Ub5_${i + 1}` : '';
            const jsonData = getMultimerData(label);

            return (
              <React.Fragment key={`${page}-${i}`}>
                {page === 'draw' ? (
                  <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', marginBottom: '16px' }}>
                    <input
                      type="text"
                      placeholder="Enter UbX_Y e.g.(Ub5_31) or A1B1B2..."
                      onKeyDown={async (e) => {
                        if (e.key === 'Enter') {
                          const value = e.target.value;
                          if (value) {
                            try {
                              // Clear previous reaction sequences and scaffolds
                              setReactionSequence(null);
                              setFigures(null);
                              setJsonOutput(null); // Clear previously rendered scaffold
                              setFormattedEdges(null); // Clear formatted edges
                              setUbxyValue(null); // Clear UbX_Y value
                              setNomenclaturePreorderABC(null); // Clear nomenclature value
                              setStritarNomenclatureWoPreorder(null); // Clear stritar nomenclature
                              setKummelstedtNomenclatureWPreorder(null); // Clear kummelstedt nomenclature
                              setInputNodes({
                                nodes: DEFAULT_NODES,
                                arrows: [],
                                mapping: {},
                                smallNodes: []
                              });
                              setScaffoldKey(Date.now());

                              const response = await fetch('/api/submit-ubxy', {
                                method: 'POST',
                                headers: { 'Content-Type': 'application/json' },
                                body: JSON.stringify({ ubxy: value }),
                              });
                              if (!response.ok) throw new Error('Failed to fetch reaction sequences');
                              const result = await response.json();
                              
                              // Store formatted edges and UbX_Y value if available
                              if (result.formatted_edges) {
                                setFormattedEdges(result.formatted_edges);
                              }
                              if (result.ubxy) {
                                setUbxyValue(result.ubxy);
                              }
                              if (result.nomenclature_preorder_ABC) {
                                setNomenclaturePreorderABC(result.nomenclature_preorder_ABC);
                              }
                              if (result.stritar_nomenclature_wo_preorder) {
                                setStritarNomenclatureWoPreorder(result.stritar_nomenclature_wo_preorder);
                              }
                              if (result.kummelstedt_nomenclature_w_preorder) {
                                setKummelstedtNomenclatureWPreorder(result.kummelstedt_nomenclature_w_preorder);
                              }
                              
                              const decodedSequence = JSON.parse(atob(result.reaction_sequences_b64));
                              setReactionSequence(decodedSequence);

                              const firstItem = decodedSequence[0];
                              const finalMultimerRaw = firstItem.ubi_his_JSON_final_multimer;
                              const finalMultimer = JSON.parse(fixQuotes(finalMultimerRaw));
                              console.log('Final Multimer:', finalMultimer);
                              setJsonOutput(finalMultimer); // Update the state for setJsonOutput

                              const simulateResult = simulateClicksFromJson(finalMultimer, DEFAULT_NODES.map(node => ({ ...node, clicks: 0 })), DEFAULT_EDGES);
                              console.log('simulateClicksFromJson output:', simulateResult);

                              // Clean the nodes and ensure all clicks are properly reset
                              let cleanedNodes;
                              cleanedNodes = simulateResult.nodes.map(({ chain_number, ...rest }) => ({ ...rest, clicks: rest.clicks || 0 }));

                              // Update inputNodes with cleaned nodes and pass all simulation data
                              setTimeout(() => {
                                setScaffoldKey(Date.now());
                                setInputNodes({
                                  nodes: cleanedNodes,
                                  arrows: simulateResult.arrows || [],
                                  mapping: simulateResult.mapping || {},
                                  smallNodes: simulateResult.smallNodes || []
                                });
                              }, 0);

                              console.log('simulateResult.nodes:', simulateResult.nodes);
                              console.log('TEST_NODES:', TEST_NODES);

                            } catch (err) {
                              alert('Failed to fetch reaction sequences: ' + err.message);
                            }
                          }
                        }
                      }}
                      style={{
                        width: '300px',
                        padding: '8px',
                        borderRadius: '6px',
                        border: '1px solid #ccc',
                        fontSize: '16px',
                        marginBottom: '12px',
                      }}
                    />
                    <div style={{
                      display: 'flex',
                      justifyContent: 'center',
                      alignItems: 'center',
                      width: panelWidth,
                      height: panelHeight,
                      margin: '0 auto',
                    }}>
                      <ClickableScaffoldPanel
                        key={scaffoldKey} // Controlled remount via state
                        initialNodes={inputNodes.nodes}
                        initialArrows={inputNodes.arrows}
                        initialMapping={inputNodes.mapping}
                        initialSmallNodes={inputNodes.smallNodes}
                        panelWidth={panelWidth}
                        panelHeight={panelHeight}
                        onSubmit={async (jsonOutput) => {
                          console.log('Linkages created:', jsonOutput);

                          // Clear previous reaction sequences and figures but preserve scaffold state
                          setReactionSequence(null);
                          setFigures(null);
                          setJsonOutput(null);
                          setFormattedEdges(null); // Clear formatted edges
                          setUbxyValue(null); // Clear UbX_Y value
                          setNomenclaturePreorderABC(null); // Clear nomenclature value
                          setStritarNomenclatureWoPreorder(null); // Clear stritar nomenclature
                          setKummelstedtNomenclatureWPreorder(null); // Clear kummelstedt nomenclature
                          // Don't reset inputNodes here to preserve the arrows and scaffold state

                          setJsonOutput(jsonOutput);
                          try {
                            const response = await fetch('/api/submit-json-output', {
                              method: 'POST',
                              headers: { 'Content-Type': 'application/json' },
                              body: JSON.stringify({ jsonOutput }),
                            });
                            if (!response.ok) {
                              console.error('Failed to fetch reaction sequences:', response.status, response.statusText);
                              throw new Error('Failed to fetch reaction sequences');
                            }
                            const result = await response.json();
                            
                            // Store formatted edges and UbX_Y value if available
                            if (result.formatted_edges) {
                              setFormattedEdges(result.formatted_edges);
                            }
                            if (result.ubxy) {
                              setUbxyValue(result.ubxy);
                            }
                            if (result.nomenclature_preorder_ABC) {
                              setNomenclaturePreorderABC(result.nomenclature_preorder_ABC);
                            }
                            if (result.stritar_nomenclature_wo_preorder) {
                              setStritarNomenclatureWoPreorder(result.stritar_nomenclature_wo_preorder);
                            }
                            if (result.kummelstedt_nomenclature_w_preorder) {
                              setKummelstedtNomenclatureWPreorder(result.kummelstedt_nomenclature_w_preorder);
                            }
                            
                            const decodedSequence = JSON.parse(atob(result.reaction_sequences_b64));
                            console.log('Decoded Sequence:', decodedSequence); // Log the decoded sequence
                            setReactionSequence(decodedSequence); // Store the fetched sequence
                          } catch (err) {
                            console.error('Error submitting jsonOutput:', err);
                          }
                        }}
                      />
                    </div>
                    {/* Display formatted edges and tree structure if available */}
                    {formattedEdges && ubxyValue && (
                      <div style={{ marginTop: '16px', width: '100%', maxWidth: '800px' }}>
                        <EdgeTreeViewer 
                          formattedEdges={formattedEdges} 
                          ubxyValue={ubxyValue} 
                          nomenclaturePreorderABC={nomenclaturePreorderABC}
                          stritarNomenclatureWoPreorder={stritarNomenclatureWoPreorder}
                          kummelstedtNomenclatureWPreorder={kummelstedtNomenclatureWPreorder}
                        />
                      </div>
                    )}
                    <div style={{ marginTop: '24px' }}> {/* Add space between ClickableScaffoldPanel and ReactionSequencesPaneled */}
                      {reactionSequence && (
                        <div style={{
                          maxHeight: '100%', // Set a fixed height for scrollability
                          overflowY: 'auto',
                          border: '1px solid #ccc',
                          borderRadius: '8px',
                          padding: '16px',
                          boxSizing: 'border-box',
                          width: '100%',
                          minWidth: '300px', // Keep the width tight
                          maxWidth: '1200px', // Keep the width tight
                          margin: '0 auto',
                        }}>
                          <ReactionSequencesPaneled reactionSequence={reactionSequence} showReactionWell={false} />
                        </div>
                      )}
                    </div>
                  </div>
                ) : (
                  <Panel
                    key={`${page}-${i}`}
                    style={{
                      width: SMALL_PANEL_WIDTH,
                      height: SMALL_PANEL_HEIGHT+35,
                      minWidth: SMALL_PANEL_WIDTH,
                      minHeight: SMALL_PANEL_HEIGHT+35,
                      padding: 0,
                      position: 'relative',
                      display: 'flex',
                      flexDirection: 'column',
                      alignItems: 'center',
                      border: isSelected ? '3px solid #1976d2' : '2px solid transparent',
                      boxShadow: isSelected ? '0 0 12px #1976d2' : undefined,
                      cursor: 'pointer',
                      transition: 'border 0.2s, box-shadow 0.2s',
                    }}
                    onClick={() => handlePanelClick(i)}
                  >
                    <div
                      style={{
                          width: SMALL_PANEL_WIDTH-5,
                          height: SMALL_PANEL_HEIGHT-2,
                          border: '1px solid #ccc',
                          borderRadius: '10px',
                          overflow: 'hidden',
                          flexShrink: 0,
                          marginBottom: '10px'
                      }}
                    >
                        <JsonToScaffold key={`${page}-${label}`} jsonData={jsonData} />
                    </div>
                    {label && (
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
                )}
              </React.Fragment>
            );
          })
        )}
      </div>
      {/* Selection tracker for tetramers and pentamers */}
      {(page === 'tetramers' || page === 'pentamers') && (
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
            <div style={{
              display: 'flex',
              flexWrap: 'wrap',
              gap: 8,
              maxWidth: '100%',
            }}>
              {Object.entries(selectionCounts).map(([label, count], idx) => (
                <span key={`selection-${label}-${idx}`} style={{
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
              ))}
            </div>
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
      {/* Render backend images if available */}
      {figures && (
        <div style={{ margin: '32px auto', maxWidth: 1200, display: 'flex', gap: 24, justifyContent: 'center', flexWrap: 'wrap' }}>
          <div style={{ width: '100%' }}>
            <button
              style={{
                marginBottom: 8,
                width: 300,
                background: '#388e3c',
                color: '#fff',
                border: 'none',
                borderRadius: 8,
                padding: '10px 0',
                fontWeight: 700,
                fontSize: 16,
                cursor: 'pointer',
                display: 'block',
                marginLeft: 'auto',
                marginRight: 'auto',
              }}
              onClick={async () => {
                const baseName = window.prompt('Enter a base file name for this group of files:', 'ubiquitin_plate');
                if (!baseName) return;
                for (const key of Object.keys(figures)) {
                  let ext = 'png';
                  let mime = 'image/png';
                  if (key.endsWith('.xlsx')) {
                    ext = 'xlsx';
                    mime = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet';
                  } else if (key.endsWith('.csv')) {
                    ext = 'csv';
                    mime = 'text/csv';
                  } else if (key.endsWith('.py')) {
                    ext = 'py';
                    mime = 'text/x-python';
                  } else if (key.endsWith('.doc') || key.endsWith('.docx')) {
                    ext = key.endsWith('.docx') ? 'docx' : 'doc';
                    mime = 'application/vnd.openxmlformats-officedocument.wordprocessingml.document';
                  } else if (key.endsWith('.txt')) {
                    ext = 'txt';
                    mime = 'text/plain';
                  } else if (key.endsWith('.json')) {
                    ext = 'json';
                    mime = 'application/json';
                  }
                  let cleanKey = key.replace(/\.(png|xlsx|csv|py|docx?|txt|json)$/i, '');
                  await new Promise(resolve => setTimeout(resolve, 200));
                  const link = document.createElement('a');
                  link.href = `data:${mime};base64,${figures[key]}`;
                  link.download = `${baseName}_${cleanKey}.${ext}`;
                  document.body.appendChild(link);
                  link.click();
                  document.body.removeChild(link);
                }
              }}
            >
              Save All Files
            </button>
            <div style={{ textAlign: 'center', color: '#388e3c', fontWeight: 600, fontSize: 15, marginBottom: 16 }}>
              The following files will also be included:<br />
              <span style={{ color: '#1976d2' }}>reagent calculations (NAME_reagent_calculations.xlsx)</span><br />
              <span style={{ color: '#1976d2' }}>opentrons synthesis (NAME_opentrons_synthesis.py)</span><br />
              <span style={{ color: '#1976d2' }}>reaction schemes (NAME_reaction_schemes.png)</span>
            </div>
          </div>
          <div>
            <div style={{ textAlign: 'center', fontWeight: 600, marginBottom: 8 }}>Enzyme + Donor Mixes</div>
            <img
              src={`data:image/png;base64,${figures.enzymes_donors_96}`}
              alt="Enzyme + Donor Plate"
              style={{ maxWidth: 350, borderRadius: 8, boxShadow: '0 2px 12px #0003' }}
            />
          </div>
          <div>
            <div style={{ textAlign: 'center', fontWeight: 600, marginBottom: 8 }}>Deprotections</div>
            <img
              src={`data:image/png;base64,${figures.deprots_96}`}
              alt="Deprotections Plate"
              style={{ maxWidth: 350, borderRadius: 8, boxShadow: '0 2px 12px #0003' }}
            />
          </div>
          <div>
            <div style={{ textAlign: 'center', fontWeight: 600, marginBottom: 8 }}>Acceptors</div>
            <img
              src={`data:image/png;base64,${figures.acceptors_96}`}
              alt="Acceptors Plate"
              style={{ maxWidth: 350, borderRadius: 8, boxShadow: '0 2px 12px #0003' }}
            />
          </div>
          <div style={{
              maxHeight: '100%',
              overflowY: 'auto',
              border: '1px solid #ccc',
              borderRadius: '8px',
              padding: '16px',
              boxSizing: 'border-box',
            }}>
            <div id="reaction-sequences">
              <ReactionSequencesPaneled reactionSequence={figures && figures["reaction_sequences.json"] ? JSON.parse(atob(figures["reaction_sequences.json"])) : null} />
            </div>          
          </div>
        </div>
      )}
    </div>
  );
};

export default ModuleDashboard;
