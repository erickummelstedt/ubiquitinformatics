import React, { useState } from 'react';
import { Bar, Scatter } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';

ChartJS.register(CategoryScale, LinearScale, BarElement, PointElement, LineElement, Title, Tooltip, Legend);

const ReactionPathStatisticsPage = () => {
  const [reactionStatsData, setReactionStatsData] = useState(null);
  const [reactionStatsLoading, setReactionStatsLoading] = useState(false);
  const [multimerSize, setMultimerSize] = useState(4);
  const [pathwayType, setPathwayType] = useState('aboc');

  // Function to fetch reaction path statistics
  const fetchReactionPathStats = async (size, type) => {
    setReactionStatsLoading(true);
    try {
      const response = await fetch('/api/reaction-path-statistics', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ multimer_size: size, pathway_type: type }),
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
        
        <div style={{ marginBottom: '24px', display: 'flex', alignItems: 'center', justifyContent: 'center', gap: '16px' }}>
          <div>
            <label style={{ marginRight: '8px', fontWeight: '600' }}>Multimer Size:</label>
            <select
              value={multimerSize}
              onChange={e => setMultimerSize(parseInt(e.target.value))}
              style={{
                padding: '8px 12px',
                borderRadius: '6px',
                border: '1px solid #ccc',
                fontSize: '16px',
                marginRight: '12px'
              }}
            >
              <option value={4}>Tetramers (4)</option>
              <option value={5}>Pentamers (5)</option>
            </select>
          </div>
          <div>
            <label style={{ marginRight: '8px', fontWeight: '600' }}>Pathway Type:</label>
            <select
              value={pathwayType}
              onChange={e => setPathwayType(e.target.value)}
              style={{
                padding: '8px 12px',
                borderRadius: '6px',
                border: '1px solid #ccc',
                fontSize: '16px',
                marginRight: '12px'
              }}
            >
              <option value="aboc">Aboc-saturated reaction pathways</option>
              <option value="all">All reaction pathways</option>
            </select>
          </div>
          <button
            onClick={() => fetchReactionPathStats(multimerSize, pathwayType)}
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
          </div>

          {/* Histogram (Bar Chart) - Number of Reactions per UbID */}
          <div style={{ 
            margin: '32px 0', 
            padding: '20px', 
            border: '1px solid #ddd', 
            borderRadius: '8px', 
            backgroundColor: 'white' 
          }}>
            <h4 style={{ margin: '0 0 16px 0', textAlign: 'center' }}>
              Histogram of Number of Reactions per UbID
            </h4>
            <div style={{ height: '400px', width: '100%', position: 'relative' }}>
              <Bar
                data={{
                  labels: Object.values(reactionStatsData.data).map(row => row.UbID),
                  datasets: [
                    {
                      label: 'Number of Reactions',
                      data: Object.values(reactionStatsData.data).map(row => row.num_of_reactions),
                      backgroundColor: '#1976d2',
                      borderColor: '#1976d2',
                      borderWidth: 1,
                    },
                  ],
                }}
                options={{
                  responsive: true,
                  maintainAspectRatio: false,
                  plugins: {
                    legend: { display: false },
                    title: { display: false },
                  },
                  scales: {
                    x: { 
                      title: { display: true, text: 'UbID' },
                      ticks: { 
                        autoSkip: false, 
                        maxRotation: 90, 
                        minRotation: 90,
                        font: { size: 10 }
                      },
                      grid: {
                        display: false
                      },
                      border: {
                        display: true,
                        width: 2,
                        color: '#333'
                      },
                      ticks: {
                        autoSkip: false,
                        maxRotation: 45,
                        minRotation: 45,
                        font: { size: 10 },
                        major: {
                          enabled: true
                        },
                        color: '#333'
                      }
                    },
                    y: { 
                      beginAtZero: true,
                      title: { display: true, text: 'Number of Reactions' },
                      grid: {
                        display: false
                      },
                      border: {
                        display: true,
                        width: 2,
                        color: '#333'
                      },
                      ticks: {
                        color: '#333'
                      }
                    },
                  },
                }}
              />
            </div>
          </div>

          {/* Scatter Plot - Number of Reactions by Linkage Type */}
          <div style={{ 
            margin: '32px 0', 
            padding: '20px', 
            border: '1px solid #ddd', 
            borderRadius: '8px', 
            backgroundColor: 'white' 
          }}>
            <h4 style={{ margin: '0 0 16px 0', textAlign: 'center' }}>
              Scatter Plot of Number of Reactions by Linkage Type
            </h4>
            <div style={{ height: '500px', width: '100%', position: 'relative' }}>
              {(() => {
                // Group data by linkage type (branching_linkage + heterotypic_linkage combination)
                const linkageGroups = {};
                Object.values(reactionStatsData.data).forEach(row => {
                  const linkageKey = `B${row.branching_linkage}_H${row.heterotypic_linkage}`;
                  if (!linkageGroups[linkageKey]) {
                    linkageGroups[linkageKey] = {
                      linkageType: linkageKey,
                      branching: row.branching_linkage,
                      heterotypic: row.heterotypic_linkage,
                      species: []
                    };
                  }
                  linkageGroups[linkageKey].species.push({
                    ubID: row.UbID,
                    reactions: row.num_of_reactions
                  });
                });

                // Sort linkage groups by branching then heterotypic (like Python code)
                const sortedLinkageTypes = Object.values(linkageGroups)
                  .sort((a, b) => {
                    if (a.branching !== b.branching) {
                      return a.branching - b.branching;
                    }
                    return a.heterotypic - b.heterotypic;
                  });

                // Create scatter plot data with tighter clustering and average lines
                const scatterData = [];
                const labels = [];
                const groupSpacing = 2; // Space between groups
                let xPosition = 0;

                sortedLinkageTypes.forEach((group, groupIndex) => {
                  labels.push(group.linkageType);
                  
                  // Constrain dots to appear only within the width of the average line
                  const lineWidth = 0.8; // Average line spans from -0.4 to +0.4
                  
                  if (group.species.length === 1) {
                    // Single point - center it exactly on the group position
                    scatterData.push({
                      x: xPosition,
                      y: group.species[0].reactions,
                      ubID: group.species[0].ubID,
                      linkageType: group.linkageType
                    });
                  } else {
                    // Multiple points - distribute them within the line width
                    const maxWidth = lineWidth;
                    const spacing = group.species.length > 1 ? maxWidth / (group.species.length - 1) : 0;
                    const startOffset = -maxWidth / 2;
                    
                    group.species.forEach((species, speciesIndex) => {
                      scatterData.push({
                        x: xPosition + startOffset + (speciesIndex * spacing),
                        y: species.reactions,
                        ubID: species.ubID,
                        linkageType: group.linkageType
                      });
                    });
                  }
                  
                  xPosition += groupSpacing; // Add spacing between groups
                });

                return (
                  <Scatter
                    data={{
                      datasets: [
                        // Create individual line datasets for each linkage group's average
                        ...sortedLinkageTypes.map((group, index) => {
                          const totalReactions = group.species.reduce((sum, species) => sum + species.reactions, 0);
                          const averageReactions = totalReactions / group.species.length;
                          return {
                            type: 'scatter',
                            label: index === 0 ? 'Average Lines' : '', // Only show legend for first one
                            data: [
                              { x: index * groupSpacing - 0.4, y: averageReactions },
                              { x: index * groupSpacing + 0.4, y: averageReactions }
                            ],
                            backgroundColor: '#ffc107',
                            borderColor: '#ffc107',
                            pointRadius: 0,
                            pointHoverRadius: 3,
                            showLine: true,
                            borderWidth: 3,
                            order: 2,
                            spanGaps: false
                          };
                        }),
                        {
                          type: 'scatter',
                          label: 'Individual Species',
                          data: scatterData,
                          backgroundColor: '#388e3c',
                          borderColor: '#2e7d32',
                          pointRadius: 6,
                          pointHoverRadius: 8,
                          showLine: false,
                          order: 1
                        }
                      ],
                    }}
                    options={{
                      responsive: true,
                      maintainAspectRatio: false,
                      layout: {
                        padding: {
                          bottom: 60
                        }
                      },
                      plugins: {
                        legend: { 
                          display: true,
                          position: 'top',
                          labels: {
                            usePointStyle: true,
                            font: { size: 12 },
                            filter: function(item) {
                              // Only show items with non-empty labels
                              return item.text && item.text.trim() !== '';
                            }
                          }
                        },
                        title: { display: false },
                        tooltip: {
                          callbacks: {
                            title: function(context) {
                              if (context[0].dataset.label === 'Individual Species') {
                                const point = scatterData[context[0].dataIndex];
                                return `${point.linkageType} - ${point.ubID}`;
                              } else {
                                // For average lines, find which group this belongs to
                                const datasetIndex = context[0].datasetIndex;
                                if (datasetIndex < sortedLinkageTypes.length) {
                                  return `${labels[datasetIndex]} - Average`;
                                }
                                return 'Average';
                              }
                            },
                            label: function(context) {
                              if (context.dataset.label === 'Individual Species') {
                                const point = scatterData[context.dataIndex];
                                return `Species: ${point.ubID}, Reactions: ${point.y}`;
                              } else {
                                return `Average Reactions: ${context.parsed.y.toFixed(1)}`;
                              }
                            }
                          }
                        }
                      },
                      scales: {
                        x: {
                          type: 'linear',
                          position: 'bottom',
                          title: { display: false },
                          min: -2,
                          max: (sortedLinkageTypes.length - 1) * groupSpacing + groupSpacing/2,
                          ticks: {
                            callback: function(value) {
                              // Multi-line headers at position -1.5 (further left)
                              if (value === -1.5) {
                                return ['Heterotypic Linkages', '', 'Branching Sites'];
                              }
                              
                              // Check if this is a center position for a linkage group
                              const mainIndex = Math.round(value / groupSpacing);
                              if (mainIndex >= 0 && mainIndex < sortedLinkageTypes.length && Math.abs(value - mainIndex * groupSpacing) < 0.1) {
                                // Return both values with extra spacing between them
                                return [
                                  sortedLinkageTypes[mainIndex].heterotypic.toString(),
                                  '', // Empty line for spacing
                                  sortedLinkageTypes[mainIndex].branching.toString()
                                ];
                              }
                              
                              return '';
                            },
                            stepSize: 0.5,
                            font: { size: 12 },
                            color: '#333',
                            maxRotation: 0,
                            minRotation: 0
                          },
                          afterBuildTicks: function(axis) {
                            // Create completely custom ticks for two-row layout with spacing
                            const customTicks = [
                              // Combined headers as multi-line label (moved further left)
                              { 
                                value: -1.5, 
                                label: ['Heterotypic Linkages', '', 'Branching Sites'] 
                              }
                            ];
                            
                            // Add ticks directly centered under each average line
                            sortedLinkageTypes.forEach((group, index) => {
                              const centerPosition = index * groupSpacing;
                              
                              // Heterotypic value at center position (first row)
                              customTicks.push({
                                value: centerPosition,
                                label: group.heterotypic.toString()
                              });
                              
                              // Branching value at center position (second row with spacing)
                              customTicks.push({
                                value: centerPosition,
                                label: group.branching.toString()
                              });
                            });
                            
                            axis.ticks = customTicks;
                          },
                          grid: {
                            display: false
                          },
                          border: {
                            display: true,
                            width: 2,
                            color: '#333'
                          }
                        },
                        y: {
                          beginAtZero: true,
                          title: { display: true, text: 'Number of Reactions' },
                          grid: {
                            display: false
                          },
                          border: {
                            display: true,
                            width: 2,
                            color: '#333'
                          },
                          ticks: {
                            color: '#333'
                          }
                        },
                      },
                    }}
                  />
                );
              })()}
            </div>
          </div>

          {/* Data Table for reference */}
          <div style={{ 
            overflowX: 'auto',
            border: '1px solid #ddd',
            borderRadius: '8px',
            backgroundColor: 'white',
            maxHeight: '50vh',
            marginTop: '32px'
          }}>
            <h4 style={{ margin: '16px', textAlign: 'center' }}>Data Table</h4>
            <table style={{ 
              width: '100%', 
              borderCollapse: 'collapse',
              fontSize: '14px'
            }}>
              <thead>
                <tr style={{ backgroundColor: '#f8f9fa' }}>
                  <th style={{ padding: '12px', textAlign: 'left', borderBottom: '2px solid #dee2e6', fontWeight: '600' }}>UbID</th>
                  <th style={{ padding: '12px', textAlign: 'left', borderBottom: '2px solid #dee2e6', fontWeight: '600' }}>Num Reactions</th>
                  <th style={{ padding: '12px', textAlign: 'left', borderBottom: '2px solid #dee2e6', fontWeight: '600' }}>UbiDAG Edges</th>
                  <th style={{ padding: '12px', textAlign: 'left', borderBottom: '2px solid #dee2e6', fontWeight: '600' }}>Heterotypic Linkage</th>
                  <th style={{ padding: '12px', textAlign: 'left', borderBottom: '2px solid #dee2e6', fontWeight: '600' }}>Branching Linkage</th>
                </tr>
              </thead>
              <tbody>
                {Object.values(reactionStatsData.data).map((row, idx) => (
                  <tr key={idx} style={{ borderBottom: '1px solid #dee2e6' }}>
                    <td style={{ padding: '12px', borderRight: '1px solid #dee2e6' }}>{row.UbID}</td>
                    <td style={{ padding: '12px', borderRight: '1px solid #dee2e6' }}>{row.num_of_reactions}</td>
                    <td style={{ padding: '12px', borderRight: '1px solid #dee2e6' }}>
                      {typeof row.ubiDAG_edges === 'string' && row.ubiDAG_edges.length > 50 
                        ? (
                          <span title={row.ubiDAG_edges}>
                            {row.ubiDAG_edges.substring(0, 50) + '...'}
                          </span>
                        )
                        : row.ubiDAG_edges}
                    </td>
                    <td style={{ padding: '12px', borderRight: '1px solid #dee2e6' }}>{row.heterotypic_linkage}</td>
                    <td style={{ padding: '12px' }}>{row.branching_linkage}</td>
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
