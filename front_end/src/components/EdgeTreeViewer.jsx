import React from 'react';

const EdgeTreeViewer = ({ 
  formattedEdges, 
  ubxyValue, 
  nomenclaturePreorder1A2, 
  nomenclaturePreorderA63B,
  graphwopreorderNomenclatureWPreorder, 
  graphwpreorderNomenclatureWPreorder,
  chemicalAllNode,
  chemicalAllNodeEricNumbering,
  outputJsonString,
  massSpecTxtFile
}) => {
  // Function to format nomenclature with subscripts for negative numbers
  const formatWithSubscripts = (text) => {
    if (!text) return text;
    
    // Replace Ub followed by numbers with subscript format
    // This will match patterns like "Ub1", "Ub2", "Ub10", etc.
    return text.replace(/Ub(\d+)/g, (match, number) => {
      // Convert each digit to its subscript Unicode equivalent
      const subscriptDigits = {
        '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
        '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
      };
      
      const subscriptNumber = number.split('').map(digit => subscriptDigits[digit]).join('');
      return `Ub${subscriptNumber}`;
    });
  };
  // Function to convert formatted edges to tree structure
  const createTreeFromEdges = (formattedEdges) => {
    if (!formattedEdges) return '';
    
    // Parse the edges string into array of connections
    const connections = formattedEdges.split(', ').map(edge => {
      const parts = edge.split(' -> ');
      return { from: parts[0], site: parts[1], to: parts[2] };
    });
    
    // Build tree structure
    const buildTree = (nodeId, depth = 0) => {
      const children = connections.filter(conn => conn.from === nodeId);
      if (children.length === 0) return '';
      
      let result = '';
      children.forEach((child, index) => {
        const isLast = index === children.length - 1;
        const indent = '          '.repeat(depth); // Increased indentation (10 spaces)
        const connector = isLast ? '└── ' : '├── ';
        result += `${indent}${connector}${child.site} → ${child.to}\n`;
        
        // Recursively build subtrees
        const subtree = buildTree(child.to, depth + 1);
        if (subtree) {
          result += subtree;
        }
      });
      
      return result;
    };
    
    // Start from node 1 (assuming it's the root)
    const treeString = '1\n' + buildTree('1');
    return treeString.trim(); // Remove any trailing newlines
  };

  if (!formattedEdges && !ubxyValue) return null;

  return (
    <div style={{
      border: '1px solid #ccc',
      borderRadius: '8px',
      padding: '16px',
      backgroundColor: '#f8f9fa',
      maxWidth: '800px',
      margin: '16px auto'
    }}>
      {/* UbX_Y Value Display */}
      {ubxyValue && (
        <div style={{ marginBottom: '16px' }}>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            color: '#333',
            fontSize: '16px',
            fontWeight: 'bold'
          }}>
            Ubiquitin Structure
          </h4>
          <div style={{
            color: '#1976d2',
            fontSize: '18px',
            fontWeight: 'bold',
            fontFamily: 'monospace',
            textAlign: 'center',
            padding: '8px',
            backgroundColor: '#e3f2fd',
            borderRadius: '4px',
            border: '1px solid #1976d2'
          }}>
            {ubxyValue}
          </div>
        </div>
      )}

      {/* Chain Linkages and Tree Structure - Display First */}
      {formattedEdges && (
        <>
          <div style={{ marginBottom: '16px' }}>
            <h4 style={{ 
              margin: '0 0 8px 0', 
              color: '#333',
              fontSize: '14px',
              fontWeight: 'bold'
            }}>
              Chain Linkages
            </h4>
            <div style={{
              color: '#555',
              fontSize: '12px',
              lineHeight: '1.3',
              fontFamily: 'monospace',
              backgroundColor: '#ffffff',
              padding: '8px',
              borderRadius: '4px',
              border: '1px solid #ddd'
            }}>
              {formattedEdges}
            </div>
          </div>
          
          <div style={{ marginBottom: '16px' }}>
            <h4 style={{ 
              margin: '0 0 8px 0', 
              color: '#333',
              fontSize: '14px',
              fontWeight: 'bold'
            }}>
              Tree Structure
            </h4>
            <div style={{
              color: '#333',
              fontSize: '12px',
              lineHeight: '1.4',
              fontFamily: 'monospace',
              textAlign: 'left',
              backgroundColor: '#f1f3f4',
              padding: '12px',
              borderRadius: '4px',
              border: '1px solid #ddd',
              whiteSpace: 'pre-wrap',
              overflowWrap: 'break-word',
              overflowX: 'auto'
            }}>
              {createTreeFromEdges(formattedEdges)}
            </div>
          </div>
        </>
      )}

      {/* Nomenclature Displays */}
      {formattedEdges && (
        <>

          {/* Ub1D Preorder Nomenclature */}
          {nomenclaturePreorderA63B && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                UbID nomenclature with preorder numbering A = 1, B = 2... (Bode/Majima/Kummelstedt)
              </h4>
              <div style={{
                color: '#1976d2',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#e3f2fd',
                borderRadius: '4px',
                border: '1px solid #1976d2'
              }}>
                {nomenclaturePreorderA63B}
              </div>
            </div>
          )}

          {/* Alternative Ub1D Preorder Nomenclature */}
          {nomenclaturePreorder1A2 && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Alternative UbID nomenclature with preorder numbering A = K63, B = K48... (Bode/Majima/Kummelstedt) 
              </h4>
              <div style={{
                color: '#7b1fa2',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#f3e5f5',
                borderRadius: '4px',
                border: '1px solid #7b1fa2'
              }}>
                {nomenclaturePreorder1A2}
              </div>
            </div>
          )}

          {/* Graph-based nomenclature with preorder numbering (Bode/Majima/Kummelstedt) */}
          {graphwpreorderNomenclatureWPreorder && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Graph-based nomenclature with preorder numbering (Bode/Majima/Kummelstedt) 
              </h4>
              <div style={{
                color: '#388e3c',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#e8f5e8',
                borderRadius: '4px',
                border: '1px solid #388e3c'
              }}>
                {formatWithSubscripts(graphwpreorderNomenclatureWPreorder)}
              </div>
            </div>
          )}

          {/* Jeff Multiple Symbols Nomenclature */}
          {chemicalAllNode && (
            <div style={{ marginBottom: '16px' }}>
              <div style={{ position: 'relative', display: 'inline-block' }}>
                <h4 style={{ 
                  margin: '0 0 8px 0', 
                  color: '#333',
                  fontSize: '14px',
                  fontWeight: 'bold',
                  cursor: 'help'
                }}
                onMouseEnter={(e) => {
                  const tooltip = e.target.nextElementSibling;
                  if (tooltip) {
                    tooltip.style.opacity = '1';
                    tooltip.style.visibility = 'visible';
                  }
                }}
                onMouseLeave={(e) => {
                  const tooltip = e.target.nextElementSibling;
                  if (tooltip) {
                    tooltip.style.opacity = '0';
                    tooltip.style.visibility = 'hidden';
                  }
                }}>
                  Chemistry-style all node labeled nomenclature (hover for info) (Bode)
                </h4>
                <div style={{
                  position: 'absolute',
                  top: '100%',
                  left: '0',
                  backgroundColor: '#333',
                  color: 'white',
                  padding: '12px',
                  borderRadius: '6px',
                  fontSize: '11px',
                  fontFamily: 'monospace',
                  lineHeight: '1.4',
                  whiteSpace: 'pre-wrap',
                  zIndex: 1000,
                  minWidth: '500px',
                  maxWidth: '600px',
                  boxShadow: '0 4px 8px rgba(0,0,0,0.2)',
                  opacity: 0,
                  visibility: 'hidden',
                  transition: 'opacity 0.3s, visibility 0.3s',
                  pointerEvents: 'none'
                }}>
                  {`Level System:
    Level 1 = A, Level 2 = B, Level 3 = C, Level 4 = D, etc.

Position Mapping:
    K63: evens with uppercase letter (e.g., B2, C2, D4)
    K48: odds with uppercase letter (e.g., B1, C1, D3)  
    K33: evens with uppercase letter + * (e.g., B*2, C*2, D*4)
    K29: odds with uppercase letter + * (e.g., B*1, C*1, D*3)
    K11: evens with lowercase letter (e.g., b2, c2, d4)
    K6:  odds with lowercase letter (e.g., b1, c1, d3)
    K27: evens with lowercase letter + * (e.g., b*2, c*2, d*4)
    M1:  odds with lowercase letter + * (e.g., b*1, c*1, d*3)

Formula: position = ((parent_letter_size + parent_number) - 1) * 2 + child_number

Where:
    parent_letter_size: Based on parent's notation type
        - 0 for uppercase without asterisk (A1, B2)
        - 2 for uppercase with asterisk (A*1, B*2)  
        - 4 for lowercase without asterisk (a1, b2)
        - 6 for lowercase with asterisk (a*1, b*2)
    parent_number: Numeric part of parent's notation
    child_number: +1 for K48/K29/K6/M1, +2 for K63/K33/K11/K27`}
                </div>
              </div>
              <div style={{
                color: '#e91e63',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#fce4ec',
                borderRadius: '4px',
                border: '1px solid #e91e63'
              }}>
                {chemicalAllNode}
              </div>
            </div>
          )}

          {/* Jeff Multiple Symbols Eric Numbering */}
          {chemicalAllNodeEricNumbering && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Jeff Multiple Symbols Eric Numbering
              </h4>
              <div style={{
                color: '#ff5722',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#fbe9e7',
                borderRadius: '4px',
                border: '1px solid #ff5722'
              }}>
                {chemicalAllNodeEricNumbering}
              </div>
            </div>
          )}

          {/* Graph-based nomenclature without preorder (Strieter/Shestoperova/Ivanov) */}
          {graphwopreorderNomenclatureWPreorder && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Graph-based nomenclature without preorder (Strieter/Shestoperova/Ivanov)
              </h4>
              <div style={{
                color: '#d32f2f',
                fontSize: '16px',
                fontWeight: 'bold',
                fontFamily: 'monospace',
                textAlign: 'center',
                padding: '8px',
                backgroundColor: '#ffebee',
                borderRadius: '4px',
                border: '1px solid #d32f2f'
              }}>
                {graphwopreorderNomenclatureWPreorder}
              </div>
            </div>
          )}

          {/* Mass spec .txt file (Strieter/Shestoperova/Ivanov) */}
          {massSpecTxtFile && (
            <div style={{ marginBottom: '16px' }}>
              <div style={{ 
                display: 'flex', 
                justifyContent: 'space-between', 
                alignItems: 'center',
                marginBottom: '8px'
              }}>
                <h4 style={{ 
                  margin: '0', 
                  color: '#333',
                  fontSize: '14px',
                  fontWeight: 'bold'
                }}>
                  Mass spec .txt file (Strieter/Shestoperova/Ivanov)
                </h4>
                <button
                  onClick={() => {
                    const blob = new Blob([massSpecTxtFile], { type: 'text/plain' });
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement('a');
                    link.href = url;
                    link.download = `${ubxyValue || 'polyubiquitin'}_mass_spec.txt`;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(url);
                  }}
                  style={{
                    backgroundColor: '#ff6f00',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    padding: '6px 12px',
                    fontSize: '12px',
                    cursor: 'pointer',
                    fontWeight: 'bold'
                  }}
                >
                  Save TXT
                </button>
              </div>
              <div style={{
                color: '#333',
                fontSize: '10px',
                lineHeight: '1.4',
                fontFamily: 'monospace',
                backgroundColor: '#fff3e0',
                padding: '12px',
                borderRadius: '4px',
                border: '1px solid #ff6f00',
                whiteSpace: 'pre-wrap',
                overflowWrap: 'break-word',
                overflowX: 'auto',
                maxHeight: '200px',
                overflowY: 'auto'
              }}>
                {massSpecTxtFile}
              </div>
            </div>
          )}
        </>
      )}

      {/* JSON Structure Display - Only show if outputJsonString is present */}
      {outputJsonString && (
        <div style={{ marginBottom: '16px' }}>
          <div style={{ 
            display: 'flex', 
            justifyContent: 'space-between', 
            alignItems: 'center',
            marginBottom: '8px'
          }}>
            <h4 style={{ 
              margin: '0', 
              color: '#333',
              fontSize: '14px',
              fontWeight: 'bold'
            }}>
              JSON Structure
            </h4>
            <button
              onClick={() => {
                const blob = new Blob([outputJsonString], { type: 'application/json' });
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.href = url;
                link.download = `${ubxyValue || 'polyubiquitin'}_structure.json`;
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                URL.revokeObjectURL(url);
              }}
              style={{
                backgroundColor: '#1976d2',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                padding: '6px 12px',
                fontSize: '12px',
                cursor: 'pointer',
                fontWeight: 'bold'
              }}
            >
              Save JSON
            </button>
          </div>
          <div style={{
            color: '#333',
            fontSize: '10px',
            lineHeight: '1.4',
            fontFamily: 'monospace',
            backgroundColor: '#f8f9fa',
            padding: '12px',
            borderRadius: '4px',
            border: '1px solid #ddd',
            whiteSpace: 'pre-wrap',
            overflowWrap: 'break-word',
            overflowX: 'auto',
            maxHeight: '300px',
            overflowY: 'auto'
          }}>
            {outputJsonString}
          </div>
        </div>
      )}
    </div>
  );
};

export default EdgeTreeViewer;
