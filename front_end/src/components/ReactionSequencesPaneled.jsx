import React, { useState, useEffect } from 'react';
import JsonToScaffold from './JsonToScaffold';


const renderBox = (item, showReactionWell) => {
    const fixQuotes = (str) => str.replace(/'/g, '"'); // Replace single quotes with double quotes

    const isPentamer = item.ubi_his_JSON_pentamer_formation; // Check if pentamer pathway exists

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
        "ubi_his_JSON_dimer_formation": item["Acceptor"],
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

    const panels = keys.map((key, index) => (
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
                            <JsonToScaffold jsonData={JSON.parse(fixQuotes(item[donorKeys[0]]))} />
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
                        <JsonToScaffold jsonData={JSON.parse(fixQuotes(item[key]))} />
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
                        <div style={{ fontSize: '24px', fontWeight: 'bold' }}>↓</div>
                    </div>
                    <div
                        style={{
                            display: 'flex',
                            flexDirection: 'column',
                            alignItems: 'center',
                            margin: '20px',
                        }}
                    >
                        <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{item["Multimer Id"]}</div>
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
                            <JsonToScaffold jsonData={JSON.parse(fixQuotes(item["ubi_his_JSON_final_multimer"]))} />
                        </div>
                    </div>
                </>
            )}
            {index > 0 && index < keys.length - 1 && (
                <>
                    {labels[key].includes('E2') && (
                        <>
                        <div
                        style={{
                            display: 'flex',
                            flexDirection: 'column',
                            alignItems: 'center',
                            margin: '0 20px',
                        }}
                        >
                            <div style={{ fontSize: '24px', fontWeight: 'bold' }}>+</div>
                        </div>
                            <div
                                style={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    alignItems: 'center',
                                    margin: '20px',
                                }}
                            >
                                <div style={{ textAlign: 'center', fontWeight: 'bold', marginBottom: '10px', whiteSpace: 'nowrap' }}>{splitReactionLabel(item[reactionLabels[index-1]], 1)}</div>
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
                                    <JsonToScaffold jsonData={JSON.parse(fixQuotes(item[donorKeys[(index/2)-1]]))} />
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
                        <div style={{ fontWeight: 'bold', marginBottom: '5px', whiteSpace: 'nowrap'  }}>{splitReactionLabel(item[reactionLabels[index-1]], 2)}</div>
                        <div style={{ fontSize: '24px', fontWeight: 'bold' }}>↓</div>
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
                            <JsonToScaffold jsonData={JSON.parse(fixQuotes(item[key]))} />
                        </div>
                    </div>
                </>
            )}
        </React.Fragment>
    ));

    return (
        <div
            style={{
                display: 'flex',
                flexDirection: 'column',
                gap: '20px',
                border: '2px solid #ccc',
                borderRadius: '10px',
                padding: '20px',
                marginBottom: '20px',
            }}
        >
            {showReactionWell && (
                <div
                    style={{
                        backgroundColor: 'rgba(255, 255, 255, 0.95)',
                        padding: '8px 12px',
                        borderRadius: '8px',
                        boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                        fontSize: '14px',
                        fontWeight: 'bold',
                        marginBottom: '10px',
                    }}
                >
                    Reaction Well #: {item["Reaction\nNumber"]}
                </div>
            )}
            <div
                style={{
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    padding: '8px 12px',
                    borderRadius: '8px',
                    boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                    fontSize: '14px',
                    fontWeight: 'bold',
                    marginBottom: '10px',
                }}
            >
                Simulation Index: {item["Simulation\nindex"]}
            </div>
            <div
                style={{
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    padding: '8px 12px',
                    borderRadius: '8px',
                    boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                    fontSize: '14px',
                    fontWeight: 'bold',
                    marginBottom: '10px',
                }}
            >
                Multimer Synthesised: {item["Multimer Id"]}
            </div>
                <div
                style={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'center',
                        margin: '20px',
                    }}
                >
                <div
                    style={{
                        width: '228px',
                        height: '148px',
                        border: '1px solid #ccc',
                        borderRadius: '10px',
                        overflow: 'hidden',
                        flexShrink: 0,
                        marginBottom: '10px'
                    }}
                >
                    <JsonToScaffold jsonData={JSON.parse(fixQuotes(item['ubi_his_JSON_final_multimer']))} />
                </div>
            </div>
            <div
                style={{
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    padding: '8px 12px',
                    borderRadius: '8px',
                    boxShadow: '0 2px 6px rgba(0,0,0,0.15)',
                    fontSize: '14px',
                    fontWeight: 'bold',
                    marginBottom: '10px',
                }}
            >
               Reaction Pathway: 
            </div>
            {panels}
        </div>
    );
};

const ReactionSequencesPaneled = ({ reactionSequence, showReactionWell = true }) => {
    const [boxes, setBoxes] = useState([]);

    useEffect(() => {
        if (reactionSequence) {
            const renderedBoxes = reactionSequence.map((item, index) => (
                <div key={item.id || index}>
                    {renderBox(item, showReactionWell)}
                </div>
            ));
            setBoxes(renderedBoxes);
        }
    }, [reactionSequence, showReactionWell]);

    return (
        <div
            className="sequences-page"
            style={{
                display: 'flex',
                flexDirection: 'row',
                gap: '20px',
                alignItems: 'flex-start',
                overflowX: 'auto',
                border: '2px solid #ccc',
                borderRadius: '10px',
                padding: '20px',
                paddingTop: '50px',
                width: '95%', // Adjust width to leave space from edges
                margin: '0 auto', // Center the content horizontally
            }}
        >
            {boxes}
        </div>
    );
};

export default ReactionSequencesPaneled;
