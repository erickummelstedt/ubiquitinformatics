import React from 'react';

const EdgeTreeViewer = ({ formattedEdges, ubxyValue, nomenclaturePreorderABC, stritarNomenclatureWoPreorder, kummelstedtNomenclatureWPreorder }) => {
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
      margin: '16px 0',
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

      {/* Formatted Edges Display */}
      {formattedEdges && (
        <>
          {/* Nomenclature Display - above Chain Linkages */}
          {nomenclaturePreorderABC && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Kummelstedt Preorder Nomenclature with ABCD... as K63, K48, K33, K29...
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
                {nomenclaturePreorderABC}
              </div>
            </div>
          )}

          {/* Stritar Nomenclature (without preorder) */}
          {stritarNomenclatureWoPreorder && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Stritar Nomenclature (without preorder)
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
                {stritarNomenclatureWoPreorder}
              </div>
            </div>
          )}

          {/* Kummelstedt Nomenclature (with preorder) */}
          {kummelstedtNomenclatureWPreorder && (
            <div style={{ marginBottom: '16px' }}>
              <h4 style={{ 
                margin: '0 0 8px 0', 
                color: '#333',
                fontSize: '14px',
                fontWeight: 'bold'
              }}>
                Kummelstedt Nomenclature (with preorder)
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
                {formatWithSubscripts(kummelstedtNomenclatureWPreorder)}
              </div>
            </div>
          )}
          
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
          
          <div>
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
    </div>
  );
};

export default EdgeTreeViewer;
