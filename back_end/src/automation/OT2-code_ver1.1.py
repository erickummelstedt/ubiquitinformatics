
from typing import List
from opentrons import protocol_api
from opentrons.protocol_api import InstrumentContext, labware
import pandas as pd
from io import StringIO
import time

metadata = {
    'protocolName': '14 tetramer synthesis',
    'author': 'Name <opentrons@example.com>',
    'description': '',
    'apiLevel': '2.13'
}


## input numbers are here ##

sample_number = 14 #write the number of ubiquitin chains to synthesize her
reaction_scale = 200
reaction_time = 960
number_of_steps = 2 #if you start from dimer and synthesize tetramer, it is 2.
isolation_time = 3 #should be 5. Only for testing.

## copy&paste of the CSV file
data = """
,1,3,5
1,SMAC_deprot,SMAC_deprot,SMAC_deprot
3,SMAC_deprot,Fake_Wash,SMAC_deprot
5,SMAC_deprot,Fake_Wash,SMAC_deprot
7,Fake_Wash,SMAC_deprot,SMAC_deprot
9,Fake_Wash,Fake_Wash,SMAC_deprot
11,Fake_Wash,SMAC_deprot,SMAC_deprot
13,Fake_Wash,SMAC_deprot,SMAC_deprot
15,Fake_Wash,SMAC_deprot,SMAC_deprot
17,SMAC_deprot,Fake_Wash,SMAC_deprot
19,Fake_Wash,SMAC_deprot,SMAC_deprot
21,Fake_Wash,SMAC_deprot,SMAC_deprot
23,SMAC_deprot,Fake_Wash,SMAC_deprot
25,SMAC_deprot,Fake_Wash,SMAC_deprot
27,SMAC_deprot,SMAC_deprot,SMAC_deprot
"""
## These are fixed variants ##
transfer_amount = 300
temperature = 40
deprotection_time = 60 #should be 60. Only for testing
speed = 200

## supporting functions
number_of_ubiquitin = int(sample_number)
number_of_lane = (number_of_ubiquitin + 8 - 1) // 8

def min(x, y):
    return x if x < y else y

def sample_to_well_positions(number_of_ubiquitin):
    rows = "ABCDEFGH"
    columns = range(1, 13)

    well_positions = []

    if number_of_ubiquitin < 1 or number_of_ubiquitin > 96:
        raise ValueError("Sample number must be between 1 and 96")

    for i in range(0, number_of_ubiquitin):
        row = rows[i % 8]
        column = columns[i // 8]
        well_positions.append(f"{row}{column}")
    
    return well_positions

well_position = sample_to_well_positions(number_of_ubiquitin)


#file_path = "2mer__to_4reaction_summary.csv"
df = pd.read_csv(StringIO(data), index_col=0)
df.shape[1] #getting column numbers
#extract_df = df.iloc[:,5::2]
#extract_df = extract_df[1::2]


## Main codes
def run(protocol: protocol_api.ProtocolContext):

    # Wash solutions from tiprack to increase the volume, using left pipette
    heater_shaker = protocol.load_module(module_name='heaterShakerModuleV1', location=1) #location is 1
    magnetic_module = protocol.load_module('magnetic module gen2', location=3)
    heating_plate = heater_shaker.load_labware('nest_96_wellplate_2ml_deep')
    magnetic_plate = magnetic_module.load_labware('nest_96_wellplate_2ml_deep')
    #wash_plate = protocol.load_labware('nest_12_reservoir_15ml', location=8)
    paradox_plate = protocol.load_labware('inhecoparadoxplate_96_tuberack_1000ul', location=10)
    tube_rack = protocol.load_labware('opentrons_6_tuberack_falcon_50ml_conical', location=8)
    #tiprack_1000 = protocol.load_labware('opentrons_96_tiprack_1000ul', location='8')
    tiprack_300_2 = protocol.load_labware('opentrons_96_tiprack_300ul', location=5)
    tiprack_300_3 = protocol.load_labware('opentrons_96_tiprack_300ul', location=6)
    tiprack_300 = protocol.load_labware('opentrons_96_tiprack_300ul', location=4)
    tiprack_300_4 = protocol.load_labware('opentrons_96_tiprack_300ul', location=9)
    left_pipette = protocol.load_instrument(
        'p300_single_gen2', mount='left', tip_racks=[tiprack_300, tiprack_300_2, tiprack_300_3, tiprack_300_4])
    right_pipette = protocol.load_instrument(
        'p300_multi_gen2', mount='right', tip_racks=[tiprack_300, tiprack_300_2, tiprack_300_3, tiprack_300_4])
    
    ## the number of the source_column[0] is the number of lanes, source_column[1] is the number of ubiquitylation (dimer = 1, trimer = 2, tetramer = 3, pentamer = 4)
    def reaction(source_column: tuple = (1, 2)):
        protocol.comment("This is the beginning of ubiquitylation reaction.")
        heater_shaker.close_labware_latch()
        right_pipette.flow_rate.aspirate = 50
        right_pipette.flow_rate.dispense = 100
        
        #add reaction mixtures to each lanes
        i = 1
        while i <= source_column[0]:
            right_pipette.pick_up_tip()
            right_pipette.transfer(
                min(reaction_scale, 300),
                paradox_plate[f'A{i + (source_column[1] - 1)*number_of_lane}'],  #source
                heating_plate[f'A{i}'],  #reaction well
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            #double if scale is more than 300 µl

            if reaction_scale > 300:
                right_pipette.transfer(
                    reaction_scale - 300,
                    paradox_plate[f'A{i + (source_column[1] - 1)*number_of_lane}'],  
                    heating_plate[f'A{i}'],
                    #This can be somewhere else. How to determine?
                    blow_out=True,
                    blowout_location='destination well',
                    new_tip='never'
                )
            right_pipette.drop_tip()
            i += 1

        heater_shaker.set_target_temperature(temperature)
        heater_shaker.set_and_wait_for_shake_speed(speed)
        protocol.delay(minutes=reaction_time)  # 16 hours reaction
        heater_shaker.deactivate_heater()
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()
        protocol.pause('reaction is done. Transfer 96 well plate to magnetic module')
    
    ## Add washing solution from tube rack
    #!!! add information on buffer (using different buffer) -> use different method
    def wash_reaction(source_column: tuple = (1, 2)):
        #removal of the solution with right pipette
        #source column[0] is number of lane, source column[1] is number of ubiquitylation, (dimer = 1, trimer = 2, tetramer = 3, pentamer = 4)
        protocol.comment("This is the beginning of wash.")
        #Removal of reaction mixture
        right_pipette.flow_rate.aspirate = 10
        magnetic_module.engage(height_from_base=4.0)
        protocol.delay(minutes=isolation_time)
        j = 1
        right_pipette.pick_up_tip()
        while j <= source_column[0]:
            right_pipette.transfer(
                reaction_scale,
                magnetic_plate[f'A{j}'],
                magnetic_plate[f'A{j + (source_column[1])*number_of_lane}'], #when preparing pentamer you have to bring this to the 
                new_tip='never'
            )
            j += 1
        right_pipette.drop_tip()
        magnetic_module.disengage()
    
        left_pipette.pick_up_tip()
        
    ###wash for 3 times
        left_pipette.flow_rate.aspirate = 150
        left_pipette.flow_rate.aspirate = 150
        for i in range(3):
            k = 0
            while k <= len(well_position) - 1:
                #left_pipette.flow_rate.aspirate = 50
                #left_pipette.pick_up_tip()
                left_pipette.transfer(
                    transfer_amount,
                    tube_rack['A1'],
                    magnetic_plate[f'{well_position[k]}'],
                    new_tip='never'
                )
                #left_pipette.mix(5, 200, magnetic_plate[f'{well_position[k]}'])
                #left_pipette.drop_tip()
                k += 1
            
            right_pipette.flow_rate.aspirate = 150
            right_pipette.flow_rate.dispense = 200

            for i in range(number_of_lane):
                right_pipette.pick_up_tip()
                right_pipette.mix(5, 200, magnetic_plate[f'A{i+1}']) #This doesn't work you have to choose A1 or A2
                right_pipette.drop_tip()

            magnetic_module.engage(height_from_base=4.0)
            protocol.delay(minutes=isolation_time)

                # wash removal using 8 channel pipettes
            l = 1
            right_pipette.pick_up_tip()
            while l <= number_of_lane:
                right_pipette.flow_rate.aspirate = 10
                right_pipette.transfer(
                    transfer_amount,
                    magnetic_plate[f'A{l}'],
                    protocol.fixed_trash['A1'],
                    new_tip='never'
                )
                l += 1
            magnetic_module.disengage()
            right_pipette.drop_tip()
        
        left_pipette.drop_tip()
        #right_pipette.drop_tip()
        protocol.comment("This is the end of wash.")
    
## deprotection method. requires source_column[0] as number_of_lane. 
## A2 = buffer, A3 = PLP solution
    def deprotection(source_column = (1)):
        protocol.comment("This is the beginning of deprotection.")
        n = 0
        left_pipette.flow_rate.aspirate = 150
        left_pipette.flow_rate.dispense = 150
        #while n <= len(well_position) - 1:
        for i in range(df.shape[0]):
            if df.iat[i,source_column] == 'SMAC_deprot':
                left_pipette.transfer(
                transfer_amount,
                tube_rack['A3'], #change later
                magnetic_plate[f'{well_position[i]}'],
                blow_out=True,
                blowout_location='destination well'
                )
            elif df.iat[i,source_column] == 'Fake_Wash':
                left_pipette.transfer(
                transfer_amount,
                tube_rack['A2'], #change later
                magnetic_plate[f'{well_position[i]}'],
                blow_out=True,
                blowout_location='destination well'
                )
            else:
                print('Error')
            #n += 1
        protocol.pause('PLP solution added. Transfer 96 well plate to heater-shaker module')
        heater_shaker.close_labware_latch()
        heater_shaker.set_target_temperature(temperature)
        heater_shaker.set_and_wait_for_shake_speed(speed)
        protocol.delay(minutes=deprotection_time)  # 60 min
        heater_shaker.deactivate_heater()
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()
        protocol.pause('Deprotection done. Transfer 96 well plate to magnetic module')
    


    def wash_deprotection(source_column: tuple = (1)):
        #removal of the solution with right pipette
        #source column[0] is number of lane

        #Removal of reaction mixture
        protocol.comment("This is the beginning of the wash of the PLP solution.")
        right_pipette.flow_rate.aspirate = 10
        right_pipette.flow_rate.dispense = 150

        magnetic_module.engage(height_from_base=4.0)
        protocol.delay(minutes=isolation_time)
        j = 1
        #right_pipette.pick_up_tip()
        while j <= source_column:
            right_pipette.transfer(
                transfer_amount,
                magnetic_plate[f'A{j}'],
                protocol.fixed_trash['A1'],
                #new_tip='never'
            )
            j += 1
        #right_pipette.drop_tip()
        magnetic_module.disengage()

        left_pipette.pick_up_tip()
        left_pipette.flow_rate.aspirate = 150
        left_pipette.flow_rate.dispense = 150
        #right_pipette.pick_up_tip()
        ### wash 3 times
        #for i in range(3):
        k = 0
        while k <= len(well_position) - 1:
            #left_pipette.flow_rate.aspirate = 50
            #left_pipette.pick_up_tip()
            left_pipette.transfer(
                transfer_amount,
                tube_rack['B1'],
                magnetic_plate[f'{well_position[k]}'],
                new_tip='never'
            )
            #left_pipette.mix(5, 200, magnetic_plate[f'{well_position[k]}'])
            #left_pipette.drop_tip()
            k += 1

        right_pipette.flow_rate.aspirate = 150
        for i in range(number_of_lane):
            right_pipette.pick_up_tip()
            right_pipette.mix(5, 200, magnetic_plate[f'A{i+1}']) #This doesn't work you have to choose A1 or A2
            right_pipette.drop_tip()
        
        magnetic_module.engage(height_from_base=4.0)
        protocol.delay(minutes=isolation_time)

            # wash removal using 8 channel pipettes
        l = 1
        right_pipette.pick_up_tip()
        while l <= number_of_lane:
            right_pipette.flow_rate.aspirate = 10
            right_pipette.transfer(
                transfer_amount,
                magnetic_plate[f'A{l}'],
                protocol.fixed_trash['A1'],
                new_tip='never'
            )
            l += 1
        magnetic_module.disengage()
        right_pipette.drop_tip()

        left_pipette.drop_tip()
        #right_pipette.drop_tip()
        protocol.pause('PLP wash done. Transfer 96 well plate to heater module')
    


    ### ending method adds buffer to the well to keep resin dry
    def ending():
        protocol.comment("This is the beginning of ending method.")
        heater_shaker.close_labware_latch()
        n = 0
        left_pipette.pick_up_tip()
        while n <= len(well_position) - 1:
            left_pipette.flow_rate.aspirate = 150
            left_pipette.transfer(
                transfer_amount,
                tube_rack['A2'],
                heating_plate[f'{well_position[n]}'],
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            n += 1
        left_pipette.drop_tip()
        heater_shaker.open_labware_latch()
        protocol.comment("This is the end of ending method.")


    current_cycle = 1
    protocol.pause('The protocol starts from Smac deprotection. Make sure to add your plate on Magnetic module')
    heater_shaker.open_labware_latch()
    while current_cycle <= number_of_steps:
        deprotection(source_column=(current_cycle-1)) #, PLP solution information
        wash_deprotection(source_column=(number_of_lane))
        reaction(source_column=(number_of_lane, current_cycle))
        wash_reaction(source_column=(number_of_lane, current_cycle))
        current_cycle += 1
    
    deprotection(source_column=(current_cycle-1))
    wash_deprotection(source_column=(number_of_lane))
    ending()
