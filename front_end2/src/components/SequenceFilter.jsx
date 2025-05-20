import React from 'react';
import Panel from './Panel';
import GameScaffoldPanel from './GameScaffoldPanel';
import Visualizer from './Visualizer';

const MOL2_URL = process.env.NODE_ENV === 'development'
  ? '/his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2'
  : 'his-GG-1ubq-1-(<K48_1ubq-2-()>).mol2';

const SequenceFilter = () => {
  return (
    <div style={{ display: 'flex', flexDirection: 'row', alignItems: 'flex-start', justifyContent: 'center', width: '100%', gap: '32px' }}>
      {/* Left/Main Panel */}
      <Panel>
        <GameScaffoldPanel />
      </Panel>
      {/* Right Hand Panel with Visualizer */}
      <Panel>
        <div style={{ display: 'flex', flexDirection: 'column', height: '100%', width: '100%', position: 'relative', overflow: 'hidden' }}>
          <div style={{ flex: 1, width: '100%', minHeight: 0, position: 'relative', overflow: 'hidden' }}>
            <Visualizer mol2Url={MOL2_URL} />
          </div>
        </div>
      </Panel>
    </div>
  );
};

export default SequenceFilter;
