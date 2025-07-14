import React, { useState, useEffect } from 'react';
import ScaffoldJsonWrapper from './ScaffoldJsonWrapper';
import reactionSequence from '/reaction_sequence.json?url';
import k48_dimer_ubiquitin from '../data/k48_dimer_ubiquitin';


const Sequences = () => {
    const [panels, setPanels] = useState([]);
    const [infoBox, setInfoBox] = useState(null);

    useEffect(() => {
        fetch(reactionSequence)
            .then((response) => response.json())
            .then((data) => {
                const fixQuotes = (str) => str.replace(/'/g, '"'); // Replace single quotes with double quotes

                const isPentamer = data.some(item => item.ubi_his_JSON_pentamer_formation); // Check if pentamer pathway exists

                const keys = [
                    "ubi_his_JSON_dimer_formation",
                    "ubi_his_JSON_dimer_deprotection",
                    "ubi_his_JSON_trimer_formation",
                    "ubi_his_JSON_trimer_deprotection",
                    "ubi_his_JSON_tetramer_formation",
                    ...(isPentamer ? ["ubi_his_JSON_tetramer_deprotection"] : []),
                    ...(isPentamer ? ["ubi_his_JSON_pentamer_formation"] : []),
                    "ubi_his_JSON_final_multimer"
                ];

                const labels = {
                    "ubi_his_JSON_dimer_formation": data[0]["Acceptor"],
                    "ubi_his_JSON_dimer_deprotection": "Deprotection:",
                    "ubi_his_JSON_trimer_formation": "E2 reaction:",
                    "ubi_his_JSON_trimer_deprotection": "Deprotection:",
                    "ubi_his_JSON_tetramer_formation": "E2 reaction:",
                    "ubi_his_JSON_tetramer_deprotection": "Deprotection:",
                    "ubi_his_JSON_pentamer_formation": "E2 reaction:",
                    "ubi_his_JSON_final_multimer": "Deprotection:"
                };

                const multimer = {
                    "ubi_his_JSON_dimer_formation": "Dimer",
                    "ubi_his_JSON_dimer_deprotection": "Dimer (deprotected)",
                    "ubi_his_JSON_trimer_formation": "Trimer",
                    "ubi_his_JSON_trimer_deprotection": "Trimer (deprotected)",
                    "ubi_his_JSON_tetramer_formation": "Tetramer",
                    "ubi_his_JSON_tetramer_deprotection": "Tetramer (deprotected)",
                    "ubi_his_JSON_pentamer_formation": "Pentamer",
                    "ubi_his_JSON_final_multimer": "Pentamer (deprotected)",
                };

                const donorKeys = [
                    "donor_JSON_trimer_formation",
                    "donor_JSON_tetramer_formation",
                    ...(isPentamer ? ["donor_JSON_pentamer_formation"] : []),
                ];

                const reactionLabels = [
                    "Dimer\ndeprotection",
                    "Trimer\nformation",
                    "Trimer\ndeprotection",
                    "Tetramer\nformation",
                    ...(isPentamer ? ["Tetramer\ndeprotection"] : []),
                    ...(isPentamer ? ["Pentamer\nformation"] : []),
                ];

                const splitReactionLabel = (label, part) => {
                    const match = label.match(/(Ubᴰ \d)\s+(.*)/);
                    if (!match) return label;
                    if (part === 1) {
                        return match[1]; // "Ubᴰ X"
                    } else if (part === 2) {
                        return match[2]; // donor enzyme
                    }
                    return label; // fallback
                };

                const mappedPanels = keys.map((key, index) => (
                    <React.Fragment key={index}>
                        {index === 0 && (
                            <div
                                style={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    alignItems: 'center',
                                    margin: '20px',
                                }}
                            >
                                {labels[key].includes('Formation') && (
                                    <div
                                        style={{
                                            width: '228px',
                                            height: '148px',
                                            border: '1px solid #ccc',
                                            borderRadius: '10px',
                                            overflow: 'hidden',
                                            flexShrink: 0,
                                            marginBottom: '10px',
                                        }}
                                    >
                                        <ScaffoldJsonWrapper jsonData={JSON.parse(fixQuotes(data[0][donorKeys[0]]))} />
                                    </div>
                                )}
                                <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{labels[key]}</div>
                                <div
                                    style={{
                                        width: '228px',
                                        height: '148px',
                                        border: '1px solid #ccc',
                                        borderRadius: '10px',
                                        overflow: 'hidden',
                                        flexShrink: 0,
                                    }}
                                >
                                    <ScaffoldJsonWrapper jsonData={JSON.parse(fixQuotes(data[0][key]))} />
                                </div>
                            </div>
                        )}
                        {index === keys.length - 1 && (
                            <>
                                <div
                                    style={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        alignItems: 'center',
                                        margin: '0 20px',
                                    }}
                                >
                                    <div style={{ fontWeight: 'bold', marginBottom: '5px', whiteSpace: 'nowrap'  }}>{labels[key]}</div>
                                    <div style={{ fontWeight: 'bold', marginBottom: '5px', whiteSpace: 'nowrap'  }}>{"Aboc"}</div>
                                    <div style={{ fontSize: '24px', fontWeight: 'bold' }}>→</div>
                                </div>
                                <div
                                    style={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        alignItems: 'center',
                                        margin: '20px',
                                    }}
                                >
                                    <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{data[0]["Multimer Id"]}</div>
                                    <div
                                        style={{
                                            width: '228px',
                                            height: '148px',
                                            border: '1px solid #ccc',
                                            borderRadius: '10px',
                                            overflow: 'hidden',
                                            flexShrink: 0,
                                        }}
                                    >
                                        <ScaffoldJsonWrapper jsonData={JSON.parse(fixQuotes(data[0][key]))} />
                                    </div>
                                </div>
                            </>
                        )}
                        {index > 0 && index < keys.length - 1 && (
                            <>
                                {labels[key].includes('E2') && (
                                    <>
                                        <div style={{ fontSize: '24px', fontWeight: 'bold' }}>+</div>
                                        <div
                                            style={{
                                                display: 'flex',
                                                flexDirection: 'column',
                                                alignItems: 'center',
                                                margin: '20px',
                                            }}
                                        >
                                            <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{splitReactionLabel(data[0][reactionLabels[index-1]], 1)}</div>
                                            <div
                                                style={{
                                                    width: '228px',
                                                    height: '148px',
                                                    border: '1px solid #ccc',
                                                    borderRadius: '10px',
                                                    overflow: 'hidden',
                                                    flexShrink: 0,
                                                }}
                                            >
                                                <ScaffoldJsonWrapper jsonData={JSON.parse(fixQuotes(data[0][donorKeys[(index/2)-1]]))} />
                                            </div>
                                        </div>
                                    </>
                                )}
                                <div
                                    style={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        alignItems: 'center',
                                        margin: '0 20px',
                                    }}
                                >
                                    <div style={{ fontWeight: 'bold', marginBottom: '5px', whiteSpace: 'nowrap'  }}>{labels[key]}</div>
                                    <div style={{ fontWeight: 'bold', marginBottom: '5px', whiteSpace: 'nowrap'  }}>{splitReactionLabel(data[0][reactionLabels[index-1]], 2)}</div>
                                    <div style={{ fontSize: '24px', fontWeight: 'bold' }}>→</div>
                                </div>
                                <div
                                    style={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        alignItems: 'center',
                                        margin: '20px',
                                    }}
                                >
                                    <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{multimer[key]}</div>
                                    <div
                                        style={{
                                            width: '228px',
                                            height: '148px',
                                            border: '1px solid #ccc',
                                            borderRadius: '10px',
                                            overflow: 'hidden',
                                            flexShrink: 0,
                                        }}
                                    >
                                        <ScaffoldJsonWrapper jsonData={JSON.parse(fixQuotes(data[0][key]))} />
                                    </div>
                                </div>
                            </>
                        )}
                    </React.Fragment>
                ));
                setPanels(mappedPanels);

                setInfoBox(
                    <div
                        style={{
                            display: 'flex',
                            flexDirection: 'row',
                            gap: '10px',
                            position: 'absolute',
                            top: '20px',
                            left: '20px',
                            zIndex: 10
                        }}
                    >
                        <div
                            style={{
                                backgroundColor: 'rgba(255, 255, 255, 0.95)',
                                padding: '8px 12px',
                                borderRadius: '8px',
                                boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                                fontSize: '14px',
                                fontWeight: 'bold'
                            }}
                        >
                            Reaction #: {data[0]["Reaction\nNumber"]}
                        </div>
                        <div
                            style={{
                                backgroundColor: 'rgba(255, 255, 255, 0.95)',
                                padding: '8px 12px',
                                borderRadius: '8px',
                                boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                                fontSize: '14px',
                                fontWeight: 'bold'
                            }}
                        >
                            Simulation Index: {data[0]["Simulation\nindex"]}
                        </div>
                    </div>
                );
            })
            .catch((error) => {
                console.error('Error fetching or parsing data:', error); // Debugging errors
            });
    }, []);

    return (
        <div
            className="sequences-page"
            style={{
                display: 'flex',
                flexDirection: 'row',
                gap: '20px',
                alignItems: 'center',
                overflowX: 'auto',
                border: '2px solid #ccc',
                borderRadius: '10px',
                padding: '20px',
                paddingTop: '50px',
                width: '95%', // Adjust width to leave space from edges
                margin: '0 auto', // Center the content horizontally
                position: 'relative'
            }}
        >
            {infoBox}
            {panels}
        </div>
    );
};

export default Sequences;
