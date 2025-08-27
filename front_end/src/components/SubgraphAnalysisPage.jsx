import React, { useState, useEffect } from 'react';

const SubgraphAnalysisPage = () => {
  const [isRunning, setIsRunning] = useState(false);
  const [csvTable, setCsvTable] = useState(null);
  const [error, setError] = useState(null);
  const [higherLevelLysineIds, setHigherLevelLysineIds] = useState(['K48', 'K63']);
  const [nLevelLysineIds, setNLevelLysineIds] = useState(['K48', 'K63']);
  const [nLevelSize, setNLevelSize] = useState(4);
  const [higherLevelSize, setHigherLevelSize] = useState(5);
  const [timingInfo, setTimingInfo] = useState(null);
  const [showTimingPopup, setShowTimingPopup] = useState(false);
  const [runningTimer, setRunningTimer] = useState(0);
  const [spinnerRotation, setSpinnerRotation] = useState(0);

  // Animation for the loading spinner
  useEffect(() => {
    let animationFrame;
    if (isRunning) {
      const animate = () => {
        setSpinnerRotation(prev => (prev + 6) % 360); // 6 degrees per frame for smooth rotation
        animationFrame = requestAnimationFrame(animate);
      };
      animationFrame = requestAnimationFrame(animate);
    }
    return () => {
      if (animationFrame) {
        cancelAnimationFrame(animationFrame);
      }
    };
  }, [isRunning]);

  const availableLysineIds = ['K6', 'K11', 'K27', 'K29', 'K33', 'K48', 'K63'];

  const formatRunningTime = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return mins > 0 ? `${mins}:${secs.toString().padStart(2, '0')}` : `${secs.toFixed(1)}s`;
  };

  const getAvailableHigherLevelSizes = () => {
    const sizes = [];
    for (let i = nLevelSize + 1; i <= 5; i++) {
      sizes.push(i);
    }
    return sizes;
  };

  const handleNLevelSizeChange = (newSize) => {
    setNLevelSize(newSize);
    // If current higher level size is not valid anymore, set it to the minimum valid size
    if (higherLevelSize <= newSize) {
      setHigherLevelSize(newSize + 1);
    }
  };

  const handleHigherLevelLysineChange = (lysineId) => {
    setHigherLevelLysineIds(prev => 
      prev.includes(lysineId) 
        ? prev.filter(id => id !== lysineId)
        : [...prev, lysineId]
    );
  };

  const handleNLevelLysineChange = (lysineId) => {
    setNLevelLysineIds(prev => {
      if (prev.includes(lysineId)) {
        return prev.filter(id => id !== lysineId);
      } else if (prev.length < 3) {
        return [...prev, lysineId];
      } else {
        setError('Maximum 3 n-level lysine IDs allowed');
        setTimeout(() => setError(null), 3000);
        return prev;
      }
    });
  };

  // Function to format edge list as individual connections
  const formatEdgeList = (edgeListStr) => {
    if (!edgeListStr) return '';
    
    console.log('formatEdgeList input:', JSON.stringify(edgeListStr)); // Debug log with JSON.stringify
    
    try {
      // Clean up the string - remove any extra whitespace and handle potential formatting issues
      let cleanStr = edgeListStr.trim();
      
      // Check if it's already in a readable format (not JSON array format)
      if (!cleanStr.startsWith('[')) {
        console.log('String does not start with [, returning as-is');
        return cleanStr;
      }
      
      // Convert single quotes to double quotes to make it valid JSON
      // This handles strings like [[1, 'K48', 2], [2, 'K48', 3]]
      cleanStr = cleanStr.replace(/'/g, '"');
      
      // Parse the edge list string - it's in format [[1, "K63", 2], [2, "K63", 3], ...]
      const edgeList = JSON.parse(cleanStr);
      if (!Array.isArray(edgeList) || edgeList.length === 0) {
        console.log('Not an array or empty array');
        return edgeListStr;
      }
      
      // Convert each edge to "from -> site -> to" format, one per line
      const formattedEdges = edgeList.map(edge => {
        if (!Array.isArray(edge) || edge.length < 3) {
          console.log('Invalid edge format:', edge);
          return edge.toString();
        }
        return `${edge[0]} -> ${edge[1]} -> ${edge[2]}`;
      });
      
      const result = formattedEdges.join('\n');
      console.log('formatEdgeList output:', JSON.stringify(result)); // Debug log
      return result;
    } catch (e) {
      console.error('Error parsing edge list:', e, 'Original string:', edgeListStr);
      // If parsing fails, return original string
      return edgeListStr;
    }
  };

  const startAnalysis = async () => {
    setIsRunning(true);
    setCsvTable(null);
    setError(null);
    setTimingInfo(null);
    setShowTimingPopup(false);
    setRunningTimer(0);
    
    // Start the timer for the running counter
    const timerInterval = setInterval(() => {
      setRunningTimer(prev => prev + 0.1);
    }, 100);
    
    // Show a "calculating estimate..." popup immediately
    setShowTimingPopup(true);
    setTimingInfo({ estimating: true });
    
    try {
      const response = await fetch('/api/analyze-subgraphs-stream', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          higher_level_size: higherLevelSize,
          n_level_size: nLevelSize,
          higher_level_lysine_ids: higherLevelLysineIds,
          n_level_lysine_ids: nLevelLysineIds
        }),
      });
      
      if (!response.ok) {
        throw new Error('Failed to start analysis');
      }

      // Handle streaming response
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let buffer = ''; // Buffer to handle partial chunks
      
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        
        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split('\n');
        
        // Keep the last line in buffer if it doesn't end with newline
        buffer = lines.pop() || '';
        
        for (const line of lines) {
          if (line.startsWith('data: ')) {
            try {
              const jsonStr = line.substring(6);
              
              // Skip empty JSON strings
              if (!jsonStr.trim()) {
                continue;
              }
              
              const data = JSON.parse(jsonStr);
              
              if (data.type === 'status') {
                // Handle status updates
                console.log('Status:', data.message);
              } else if (data.type === 'progress_update') {
                // Handle progress updates
                console.log(`Progress: ${data.data.current}/${data.data.total} - ${data.data.message}`);
              } else if (data.type === 'timing_update') {
                // Handle timing updates - this is the key fix!
                setTimingInfo(data.data);
                setShowTimingPopup(true);
                // Don't auto-hide for timing updates - let user manually close
              } else if (data.type === 'final_results') {
                // Handle final results
                const result = data.data;
                
                // Handle final timing information
                if (result.timing_analysis) {
                  setTimingInfo(result.timing_analysis);
                  setShowTimingPopup(true);
                  // Auto-hide the popup after 10 seconds for final results
                  setTimeout(() => setShowTimingPopup(false), 10000);
                } else {
                  // If no timing info, hide the estimating popup
                  setShowTimingPopup(false);
                }
                
                // Decode base64 CSV and parse into table
                if (result.csv_b64) {
                  const csvString = atob(result.csv_b64);
                  const rows = csvString.trim().split('\n').map(row => {
                    // Proper CSV parsing to handle commas inside quoted strings
                    const cells = [];
                    let current = '';
                    let inQuotes = false;
                    
                    for (let i = 0; i < row.length; i++) {
                      const char = row[i];
                      if (char === '"') {
                        inQuotes = !inQuotes;
                        current += char;
                      } else if (char === ',' && !inQuotes) {
                        cells.push(current);
                        current = '';
                      } else {
                        current += char;
                      }
                    }
                    cells.push(current); // Add the last cell
                    return cells;
                  });
                  setCsvTable(rows);
                } else {
                  setError('No CSV data received from API');
                }
              } else if (data.type === 'error') {
                setError(data.message);
                setShowTimingPopup(false);
              }
            } catch (parseError) {
              console.error('Error parsing stream data:', parseError, 'Line:', line);
              // If we get consistent JSON parse errors, fall back to non-streaming
              if (parseError.message.includes('Unterminated string') || parseError.message.includes('JSON Parse error')) {
                console.log('Falling back to non-streaming endpoint due to JSON parse errors');
                setError('Streaming failed, falling back to standard analysis...');
                return startAnalysisNonStreaming();
              }
            }
          }
        }
      }
    } catch (err) {
      setError(err.message);
      setShowTimingPopup(false); // Hide popup on error
    } finally {
      clearInterval(timerInterval);
      setIsRunning(false);
    }
  };

  const startAnalysisNonStreaming = async () => {
    // Start the timer for the running counter
    const timerInterval = setInterval(() => {
      setRunningTimer(prev => prev + 0.1);
    }, 100);
    
    try {
      const response = await fetch('/api/analyze-subgraphs', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          higher_level_size: higherLevelSize,
          n_level_size: nLevelSize,
          higher_level_lysine_ids: higherLevelLysineIds,
          n_level_lysine_ids: nLevelLysineIds
        }),
      });
      
      if (!response.ok) {
        throw new Error('Failed to start analysis');
      }
      
      const result = await response.json();
      
      // Handle timing information if available
      if (result.timing_analysis) {
        setTimingInfo(result.timing_analysis);
        setShowTimingPopup(true);
        // Auto-hide the popup after 10 seconds since this shows the actual estimate
        setTimeout(() => setShowTimingPopup(false), 10000);
      } else {
        // If no timing info, hide the estimating popup
        setShowTimingPopup(false);
      }
      
      // Decode base64 CSV and parse into table
      if (result.csv_b64) {
        const csvString = atob(result.csv_b64);
        const rows = csvString.trim().split('\n').map(row => {
          // Proper CSV parsing to handle commas inside quoted strings
          const cells = [];
          let current = '';
          let inQuotes = false;
          
          for (let i = 0; i < row.length; i++) {
            const char = row[i];
            if (char === '"') {
              inQuotes = !inQuotes;
              current += char;
            } else if (char === ',' && !inQuotes) {
              cells.push(current);
              current = '';
            } else {
              current += char;
            }
          }
          cells.push(current); // Add the last cell
          return cells;
        });
        setCsvTable(rows);
      } else {
        setError('No CSV data received from API');
      }
    } catch (err) {
      setError(err.message);
      setShowTimingPopup(false); // Hide popup on error
    } finally {
      clearInterval(timerInterval);
    }
  };

  return (
    <div style={{ padding: '32px', fontFamily: 'Arial, sans-serif' }}>
      <h2>Ubiquitin Chain Isomorphisms</h2>
      <div style={{ marginBottom: '24px', padding: '16px', backgroundColor: '#f0f8ff', border: '1px solid #b0d4f1', borderRadius: '8px' }}>
        <h3 style={{ margin: '0 0 12px 0', color: '#2c5282', fontWeight: 'bold' }}>Subgraph Containment Analysis Table</h3>
        <p style={{ margin: 0, color: '#4a5568', lineHeight: '1.5' }}>
          This table displays a subgraph containment matrix that quantifies how many times smaller ubiquitin multimers (subgraphs) are found as isomorphic structures within larger ubiquitin multimers (supergraphs). Each subgraph corresponds to a structural epitope.
        </p>
      </div>
      
      {/* Configuration Section */}
      <div style={{ marginBottom: '32px', padding: '20px', border: '1px solid #ddd', borderRadius: '8px', backgroundColor: '#f9f9f9' }}>
        <h3 style={{ marginTop: 0, marginBottom: '20px' }}>Analysis Configuration</h3>
        
        {/* Size Selection */}
        <div style={{ marginBottom: '20px' }}>
          <h4 style={{ marginBottom: '10px' }}>Multimer Sizes:</h4>
          <div style={{ display: 'flex', gap: '20px', alignItems: 'center', marginBottom: '10px' }}>
            <div>
              <label style={{ marginRight: '10px', fontWeight: 'bold' }}>N-Level Size (2-4):</label>
              <select 
                value={nLevelSize} 
                onChange={(e) => handleNLevelSizeChange(parseInt(e.target.value))}
                style={{ padding: '5px', fontSize: '14px', borderRadius: '4px', border: '1px solid #ccc' }}
              >
                <option value={2}>2 (Dimers)</option>
                <option value={3}>3 (Trimers)</option>
                <option value={4}>4 (Tetramers)</option>
              </select>
            </div>
            <div>
              <label style={{ marginRight: '10px', fontWeight: 'bold' }}>Higher Level Size:</label>
              <select 
                value={higherLevelSize} 
                onChange={(e) => setHigherLevelSize(parseInt(e.target.value))}
                style={{ padding: '5px', fontSize: '14px', borderRadius: '4px', border: '1px solid #ccc' }}
              >
                {getAvailableHigherLevelSizes().map(size => (
                  <option key={size} value={size}>
                    {size} ({size === 2 ? 'Dimers' : size === 3 ? 'Trimers' : size === 4 ? 'Tetramers' : 'Pentamers'})
                  </option>
                ))}
              </select>
            </div>
          </div>
          <div style={{ fontSize: '12px', color: '#666', fontStyle: 'italic' }}>
            Higher level size must be greater than n-level size
          </div>
        </div>
        
        {/* N-Level Lysine IDs */}
        <div style={{ marginBottom: '20px' }}>
          <h4 style={{ marginBottom: '10px' }}>N-Level Lysine IDs ({nLevelSize === 2 ? 'Dimers' : nLevelSize === 3 ? 'Trimers' : 'Tetramers'}) - Max 3:</h4>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '10px' }}>
            {availableLysineIds.map(lysineId => (
              <label key={`n-${lysineId}`} style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                <input
                  type="checkbox"
                  checked={nLevelLysineIds.includes(lysineId)}
                  onChange={() => handleNLevelLysineChange(lysineId)}
                  disabled={!nLevelLysineIds.includes(lysineId) && nLevelLysineIds.length >= 3}
                  style={{ marginRight: '5px' }}
                />
                {lysineId}
              </label>
            ))}
          </div>
          <div style={{ marginTop: '5px', fontSize: '12px', color: '#666' }}>
            Selected: {nLevelLysineIds.join(', ') || 'None'} ({nLevelLysineIds.length}/3)
          </div>
        </div>
        
        {/* Higher Level Lysine IDs */}
        <div style={{ marginBottom: '20px' }}>
          <h4 style={{ marginBottom: '10px' }}>Higher Level Lysine IDs ({higherLevelSize === 2 ? 'Dimers' : higherLevelSize === 3 ? 'Trimers' : higherLevelSize === 4 ? 'Tetramers' : 'Pentamers'}):</h4>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '10px' }}>
            {availableLysineIds.map(lysineId => (
              <label key={`higher-${lysineId}`} style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                <input
                  type="checkbox"
                  checked={higherLevelLysineIds.includes(lysineId)}
                  onChange={() => handleHigherLevelLysineChange(lysineId)}
                  style={{ marginRight: '5px' }}
                />
                {lysineId}
              </label>
            ))}
          </div>
          <div style={{ marginTop: '5px', fontSize: '12px', color: '#666' }}>
            Selected: {higherLevelLysineIds.join(', ') || 'None'}
          </div>
        </div>
      </div>

      <button
        onClick={startAnalysis}
        disabled={isRunning || higherLevelLysineIds.length === 0 || nLevelLysineIds.length === 0}
        style={{
          padding: '12px 24px',
          backgroundColor: (isRunning || higherLevelLysineIds.length === 0 || nLevelLysineIds.length === 0) ? '#94a3b8' : '#059669',
          color: 'white',
          border: 'none',
          borderRadius: '6px',
          cursor: (isRunning || higherLevelLysineIds.length === 0 || nLevelLysineIds.length === 0) ? 'not-allowed' : 'pointer',
          fontSize: '16px',
          fontWeight: 'bold',
          marginBottom: '24px',
          display: 'flex',
          alignItems: 'center',
          gap: '8px',
          transition: 'background-color 0.2s ease'
        }}
      >
        {isRunning ? (
          <>
            <div style={{
              width: '16px',
              height: '16px',
              border: '2px solid rgba(255, 255, 255, 0.3)',
              borderTop: '2px solid white',
              borderRadius: '50%',
              transform: `rotate(${spinnerRotation}deg)`
            }}></div>
            {`Running... ${formatRunningTime(runningTimer)}`}
          </>
        ) : (
          'Start Subgraph Analysis'
        )}
      </button>
      {error && (
        <div style={{ color: 'red', marginBottom: '16px' }}>Error: {error}</div>
      )}
      
      {/* Timing Information Popup */}
      {showTimingPopup && timingInfo && (
        <div style={{
          position: 'fixed',
          top: '20px',
          right: '20px',
          backgroundColor: timingInfo.estimating ? '#fff3cd' : '#e8f5e8',
          border: `2px solid ${timingInfo.estimating ? '#ffeaa7' : '#4caf50'}`,
          borderRadius: '8px',
          padding: '16px',
          boxShadow: '0 4px 12px rgba(0,0,0,0.15)',
          zIndex: 1000,
          minWidth: '300px',
          maxWidth: '400px'
        }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '12px' }}>
            <h4 style={{ margin: 0, color: timingInfo.estimating ? '#856404' : '#2e7d32', fontWeight: 'bold' }}>
              {timingInfo.estimating ? 'Calculating Time Estimate...' : 'Analysis Timing Estimate'}
            </h4>
            <button
              onClick={() => setShowTimingPopup(false)}
              style={{
                background: 'none',
                border: 'none',
                fontSize: '18px',
                cursor: 'pointer',
                color: '#666',
                padding: '0',
                marginLeft: '10px'
              }}
            >
              ×
            </button>
          </div>
          <div style={{ fontSize: '14px', lineHeight: '1.4', color: '#333' }}>
            {timingInfo.estimating ? (
              <div>
                <div style={{ marginBottom: '8px' }}>
                  Analysis is running... Time estimate will be available after the first 10 iterations.
                </div>
                <div style={{ fontSize: '12px', color: '#666', fontStyle: 'italic' }}>
                  This helps predict total completion time based on initial performance sampling.
                </div>
              </div>
            ) : (
              <>
                <div style={{ marginBottom: '8px' }}>
                  <strong>Estimated Total Time:</strong> {timingInfo.estimated_total_seconds ? 
                    `${timingInfo.estimated_total_seconds.toFixed(1)} seconds` : 
                    'Calculating...'}
                </div>
                {timingInfo.estimated_remaining_seconds && (
                  <div style={{ marginBottom: '8px' }}>
                    <strong>Estimated Remaining:</strong> {timingInfo.estimated_remaining_seconds.toFixed(1)} seconds
                  </div>
                )}
                <div style={{ marginBottom: '8px' }}>
                  <strong>Avg Time per Iteration:</strong> {timingInfo.avg_time_per_iteration ? 
                    `${timingInfo.avg_time_per_iteration.toFixed(2)} seconds` : 
                    'Calculating...'}
                </div>
                {timingInfo.completed_iterations && (
                  <div style={{ fontSize: '12px', color: '#666', fontStyle: 'italic' }}>
                    Based on analysis of first {timingInfo.completed_iterations} iterations
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      )}
      
      {csvTable && (
        <div style={{ marginTop: '16px', width: '100%', maxWidth: '100vw' }}>
          <h4>Results Table:</h4>
          <div style={{ 
            overflowX: 'auto', 
            overflowY: 'auto', 
            background: '#f4f4f4', 
            padding: '12px', 
            borderRadius: '4px',
            maxHeight: '70vh',
            width: '100%',
            border: '1px solid #ddd'
          }}>
            <table style={{ 
              borderCollapse: 'collapse', 
              width: 'max-content',
              minWidth: '100%'
            }}>
              <thead style={{ position: 'sticky', top: 0, zIndex: 10 }}>
                <tr>
                  {csvTable[0].map((cell, idx) => {
                    // Extract content between double quotes, if present
                    const match = cell.match(/"([^"]*)"/);
                    const displayHeader = match ? match[1] : cell;
                    const formattedHeader = formatEdgeList(displayHeader);
                    return (
                      <th key={idx} style={{ 
                        border: '1px solid #ccc', 
                        padding: '8px', 
                        background: '#e0e0e0', 
                        fontWeight: 'bold', 
                        fontSize: '11px', 
                        whiteSpace: 'pre-wrap', 
                        fontFamily: 'monospace', 
                        minWidth: '100px', 
                        maxWidth: '150px', 
                        wordWrap: 'break-word', 
                        textAlign: 'center',
                        position: 'sticky',
                        top: 0
                      }}>
                        {formattedHeader}
                      </th>
                    );
                  })}
                </tr>
              </thead>
              <tbody>
                {csvTable.slice(1).map((row, rIdx) => (
                  <tr key={rIdx}>
                    {row.map((cell, cIdx) => {
                      if (cIdx === 0) {
                        // First column is row header - format as edge list
                        const match = cell.match(/"([^"]*)"/);
                        const displayHeader = match ? match[1] : cell;
                        const formattedHeader = formatEdgeList(displayHeader);
                        return (
                          <th key={cIdx} style={{ 
                            border: '1px solid #ccc', 
                            padding: '6px', 
                            fontSize: '10px', 
                            background: '#e0e0e0', 
                            fontWeight: 'bold', 
                            whiteSpace: 'pre-wrap', 
                            fontFamily: 'monospace', 
                            textAlign: 'center',
                            minWidth: '100px',
                            maxWidth: '150px',
                            wordWrap: 'break-word',
                            position: 'sticky',
                            left: 0,
                            zIndex: 5
                          }}>
                            {formattedHeader}
                          </th>
                        );
                      } else {
                        // Regular data cell - check if value is non-zero to make it bold
                        const cellValue = cell.trim();
                        const isNonZero = cellValue !== '0' && cellValue !== '0.0' && cellValue !== '' && cellValue !== '0.00';
                        return (
                          <td key={cIdx} style={{ 
                            border: '1px solid #ccc', 
                            padding: '4px', 
                            fontSize: '12px', 
                            background: rIdx % 2 === 0 ? '#fff' : '#f9f9f9',
                            textAlign: 'center',
                            fontWeight: isNonZero ? 'bold' : 'normal',
                            minWidth: '60px',
                            maxWidth: '100px'
                          }}>
                            {cell}
                          </td>
                        );
                      }
                    })}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
};

export default SubgraphAnalysisPage;
