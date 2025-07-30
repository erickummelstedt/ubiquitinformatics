import React, { useState } from 'react';

const ReactionPathStatisticsPage = () => {
  const [reactionStatsData, setReactionStatsData] = useState(null);
  const [reactionStatsLoading, setReactionStatsLoading] = useState(false);

  // Function to fetch reaction path statistics
  const fetchReactionPathStats = async (multimerSize) => {
    setReactionStatsLoading(true);
    try {
      const response = await fetch('/api/reaction-path-statistics', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ multimer_size: multimerSize }),
      });
      
      if (!response.ok) {
        // Get the error details from the response
        const errorData = await response.text();
        console.error('Response error:', errorData);
        throw new Error(`HTTP ${response.status}: ${errorData}`);
      }
      
      const result = await response.json();
      console.log('Full API response:', result);
      
      // Check if the response has an error status
      if (result.status === 'error') {
        throw new Error(result.message || 'Unknown server error');
      }
      
      // Alert the json_with_reaction_information data
      alert('Reaction Path Statistics Data:\n' + JSON.stringify(result.data, null, 2));
      
      setReactionStatsData(result);
    } catch (err) {
      console.error('Full error details:', err);
      alert('Failed to fetch reaction path statistics: ' + err.message);
    } finally {
      setReactionStatsLoading(false);
    }
  };

  return (
    <div>
      <div style={{ textAlign: 'center', marginBottom: '24px' }}>
        <h2 style={{ margin: '0 0 16px 0', color: '#333' }}>Reaction Path Statistics</h2>
        <p style={{ margin: '0 0 24px 0', color: '#666' }}>
          Analyze reaction pathways and linkage patterns for multimers
        </p>
        
        <div style={{ marginBottom: '24px' }}>
          <label style={{ marginRight: '12px', fontWeight: '600' }}>Multimer Size:</label>
          <select 
            defaultValue="5"
            onChange={(e) => fetchReactionPathStats(parseInt(e.target.value))}
            style={{ 
              padding: '8px 12px', 
              borderRadius: '6px', 
              border: '1px solid #ccc',
              fontSize: '16px',
              marginRight: '12px'
            }}
          >
            <option value="4">Tetramers (4)</option>
            <option value="5">Pentamers (5)</option>
          </select>
          
          <button
            onClick={() => fetchReactionPathStats(5)}
            disabled={reactionStatsLoading}
            style={{
              padding: '8px 16px',
              backgroundColor: reactionStatsLoading ? '#ccc' : '#1976d2',
              color: 'white',
              border: 'none',
              borderRadius: '6px',
              cursor: reactionStatsLoading ? 'not-allowed' : 'pointer',
              fontSize: '16px'
            }}
          >
            {reactionStatsLoading ? 'Loading...' : 'Analyze'}
          </button>
        </div>
      </div>

      {reactionStatsData && (
        <div>
          <div style={{ 
            backgroundColor: '#f5f5f5', 
            padding: '16px', 
            borderRadius: '8px', 
            marginBottom: '24px',
            textAlign: 'center'
          }}>
            <h3 style={{ margin: '0 0 8px 0' }}>Analysis Results</h3>
            <p style={{ margin: '0', color: '#666' }}>
              {reactionStatsData.message}
            </p>
            {reactionStatsData.csv_b64 && (
              <button
                onClick={() => {
                  const link = document.createElement('a');
                  link.href = `data:text/csv;base64,${reactionStatsData.csv_b64}`;
                  link.download = `reaction_path_statistics_size_${reactionStatsData.multimer_size}.csv`;
                  document.body.appendChild(link);
                  link.click();
                  document.body.removeChild(link);
                }}
                style={{
                  marginTop: '12px',
                  padding: '8px 16px',
                  backgroundColor: '#388e3c',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '14px'
                }}
              >
                Download CSV
              </button>
            )}
          </div>

          <div style={{ 
            overflowX: 'auto',
            border: '1px solid #ddd',
            borderRadius: '8px',
            backgroundColor: 'white',
            maxHeight: '70vh'
          }}>
            <table style={{ 
              width: '100%', 
              borderCollapse: 'collapse',
              fontSize: '14px'
            }}>
              <thead>
                <tr style={{ backgroundColor: '#f8f9fa' }}>
                  {reactionStatsData.columns && reactionStatsData.columns.map(col => (
                    <th key={col} style={{ 
                      padding: '12px', 
                      textAlign: 'left', 
                      borderBottom: '2px solid #dee2e6',
                      fontWeight: '600',
                      position: 'sticky',
                      top: 0,
                      backgroundColor: '#f8f9fa',
                      zIndex: 10
                    }}>
                      {col}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {reactionStatsData.data && reactionStatsData.data.map((row, idx) => (
                  <tr key={idx} style={{ 
                    borderBottom: '1px solid #dee2e6'
                  }}>
                    {reactionStatsData.columns.map(col => (
                      <td key={col} style={{ 
                        padding: '12px',
                        borderRight: '1px solid #dee2e6',
                        verticalAlign: 'top'
                      }}>
                        {typeof row[col] === 'string' && row[col].length > 50 
                          ? (
                            <span title={row[col]}>
                              {row[col].substring(0, 50) + '...'}
                            </span>
                          )
                          : row[col]}
                      </td>
                    ))}
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

export default ReactionPathStatisticsPage;
