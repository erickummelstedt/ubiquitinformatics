#transfer everything on the magnetic module
# personal memo
# update from this code
# default condition is changed from tetramer to pentamer - done
# reaction time is 630 min - done
# centrifuge was added after opening the heater shaker - done
# centrifuge was added after deprotection - done
# switch the last equilibration to simply adding wash buffer
# change the height of pipette when dispensing


from typing import List
from opentrons import protocol_api
from opentrons.protocol_api import InstrumentContext, labware
from io import StringIO
import time

metadata = {
    "ctxName": "OTflex_dimer formation",
    "author": "Sohei Majima",
    "description": "",
}

requirements = {"robotType": "Flex", "apiLevel": "2.18"}

##############################


## copy&paste of the CSV file

data = [['SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot'],
 ['Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash',
  'Fake_Wash'],
 ['SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot',
  'SMAC_deprot']]

###################

isolation_time = 1.5 #should be 5. Only for testing.
transfer_amount = 200
temperature = 40
speed = 700


def add_parameters(parameters):
    parameters.add_int(
        variable_name="num_sample",
        display_name="number_of_sample",
        description="How many wells do you need?",
        default=16,
        minimum=1,
        maximum=16
        )
    
    parameters.add_int(
        variable_name="ub_len",
        display_name="ub_length",
        description="Are you synthesizing tetramer?",
        default=5,
        choices=[
            {"display_name": "tetramer", "value": 4},
            {"display_name": "pentamer", "value": 5},
            {"display_name": "dimer", "value": 2},
        ],
    )

    parameters.add_int(
        variable_name="scale",
        display_name="scale",
        description="What's the scale of the reaction in µl?",
        default=200,
        minimum=100,
        maximum=600
        )
    
    parameters.add_int(
        variable_name="acceptor",
        display_name="acceptor",
        description="What is your acceptor?",
        default=2,
        choices=[
            {"display_name": "dimer", "value": 2},
            {"display_name": "monomer", "value": 1},
        ],
    )

    parameters.add_int(
        variable_name="test",
        display_name="test",
        description="Is this a test run (short reaction time)?",
        default=0,
        choices=[
            {"display_name": "yes", "value": 1},
            {"display_name": "no", "value": 0},
        ],
    )

    parameters.add_int(
        variable_name="shaker",
        display_name="shaker speed",
        description="What's the speed of the shaker for wash?",
        default=800,
        minimum=100,
        maximum=1000
        )

    

##########################


def run(ctx):

    global num_sample
    global ub_len
    global scale
    global acceptor
    global test
    global shaker

    num_sample = ctx.params.num_sample
    ub_len = ctx.params.ub_len
    scale = ctx.params.scale   
    acceptor = ctx.params.acceptor
    test = ctx.params.test
    shaker = ctx.params.shaker

    if test == 1:
        reaction_time = 0.5
        deprotection_time = 0.5
    else:
        reaction_time = 600 #stick to it if you want to run it with in 3 days
        deprotection_time = 120


    # load labware
    heater_shaker = ctx.load_module('heaterShakerModuleV1', 'C1')
    hs_adapter = heater_shaker.load_adapter('opentrons_96_deep_well_adapter')
    magnetic_module = ctx.load_module("magneticBlockV1", 'D1')
    mm_adapter = magnetic_module.load_adapter('opentrons_96_deep_well_adapter')
    reaction_plate = mm_adapter.load_labware('nest_96_wellplate_2ml_deep')
    tube_rack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', "D2")
    tiprack_200 = ctx.load_labware("opentrons_flex_96_tiprack_200ul", "B3")
    tiprack_200_2 = ctx.load_labware("opentrons_flex_96_tiprack_200ul", "C3")
    tiprack_200_3 = ctx.load_labware("opentrons_flex_96_tiprack_200ul", "D3")
    tiprack_200_4 = ctx.load_labware("opentrons_flex_96_tiprack_200ul", "A2")
    tiprack_200_5 = ctx.load_labware("opentrons_flex_96_tiprack_200ul", "B1")
    paradox_plate = ctx.load_labware('inhecoparadoxplate_96_tuberack_1000ul', "A1")

    trash = ctx.load_trash_bin('A3')

    left_pipette = ctx.load_instrument(
        'flex_1channel_1000', 'left', tip_racks=[tiprack_200, tiprack_200_2, tiprack_200_3, tiprack_200_4, tiprack_200_5])
    right_pipette = ctx.load_instrument(
        'flex_8channel_1000', 'right', tip_racks=[tiprack_200, tiprack_200_2, tiprack_200_3, tiprack_200_4, tiprack_200_5])

    default_rate = 300
    left_pipette.flow_rate.aspirate = default_rate
    right_pipette.flow_rate.dispense = default_rate

    
    wash_stock = ctx.load_labware('nest_12_reservoir_15ml', 'B2', 'wash buffer')
    depro_wash_stock = ctx.load_labware('nest_12_reservoir_15ml', 'C2', 'deprotection wash buffer')

    ####liquid
    wash_def = ctx.define_liquid(name="WASH", description="Wash Buffer", display_color="#9ACECB")  

    depro_wash_def = ctx.define_liquid(name="DEPROTECTION_WASH", description="Deprotection Wash Buffer", display_color="#808080")  

    wash_stock["A1"].load_liquid(liquid=wash_def, volume=15000) 
    depro_wash_stock["A1"].load_liquid(liquid=depro_wash_def, volume=15000) 

    ## supporting functions

    number_of_ubiquitin = int(num_sample)
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

    #########
    def reaction(source_column: tuple = (1, 2)):
        ctx.comment("This is the beginning of ubiquitylation reaction.")
        right_pipette.flow_rate.aspirate = 50
        right_pipette.flow_rate.dispense = 10
        
        #add reaction mixtures to each lanes
        i = 1
        while i <= source_column[0]:
            right_pipette.pick_up_tip()
            right_pipette.transfer(
                min(scale, 200),
                paradox_plate[f'A{i + (source_column[1] - 1)*number_of_lane}'],  #source
                reaction_plate[f'A{i}'].bottom(3),  #reaction well
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            #double if scale is more than 300 µl

            if scale > 200:
                right_pipette.transfer(
                    scale - 200,
                    paradox_plate[f'A{i + (source_column[1] - 1)*number_of_lane}'],  
                    reaction_plate[f'A{i}'].bottom(3),
                    #This can be somewhere else. How to determine?
                    blow_out=True,
                    blowout_location='destination well',
                    new_tip='never'
                )
            right_pipette.drop_tip()
            i += 1
        
        # adding lid
        ctx.move_labware(labware = reaction_plate,
                              new_location=hs_adapter,
                              use_gripper=True
                              )


        heater_shaker.close_labware_latch()
        ctx.pause('Seal the 96 well plate')

        heater_shaker.set_target_temperature(temperature)
        heater_shaker.set_and_wait_for_shake_speed(speed)
        ctx.delay(minutes=reaction_time)  # 16 hours reaction
        heater_shaker.deactivate_heater()
        heater_shaker.deactivate_shaker()

    
    def wash_reaction(source_column: tuple = (1, 2)):
        #removal of the solution with right pipette
        #source column[0] is number of lane, source column[1] is number of ubiquitylation, (dimer = 1, trimer = 2, tetramer = 3, pentamer = 4)
        ctx.comment("This is the beginning of wash.")
        #Removal of reaction mixture
        right_pipette.flow_rate.aspirate = 10
        right_pipette.flow_rate.dispense = 200
        heater_shaker.open_labware_latch()
        # removing lid
        ctx.pause('Remove the seal of the 96 well plate and centrifuge the plate')
        #moving from heater shaker to magnetic module
        ctx.move_labware(labware = reaction_plate,
                              new_location=mm_adapter,
                              use_gripper=True
                              )
        ctx.delay(minutes=isolation_time)
        j = 1
        right_pipette.pick_up_tip()
        while j <= source_column[0]:
            right_pipette.transfer(
                scale,
                reaction_plate[f'A{j}'],
                reaction_plate[f'A{j + (source_column[1])*number_of_lane}'], 
                new_tip='never'
            )
            j += 1
        right_pipette.drop_tip()
        
    ###wash for 3 times
        for i in range(2):
            right_pipette.flow_rate.aspirate = 200
            right_pipette.flow_rate.dispense = 10
            right_pipette.pick_up_tip()
            k = 0
            while k <= number_of_lane - 1:
                #left_pipette.flow_rate.aspirate = 50
                #left_pipette.pick_up_tip()
                right_pipette.transfer(
                    transfer_amount,
                    wash_stock['A1'],
                    reaction_plate[f'A{k+1}'].bottom(3),
                    new_tip='never'
                )
                #left_pipette.mix(5, 200, magnetic_plate[f'{well_position[k]}'])
                #left_pipette.drop_tip()
                k += 1
            
            right_pipette.drop_tip()

            ctx.move_labware(labware = reaction_plate,
                              new_location=hs_adapter,
                              use_gripper=True
                              )
            heater_shaker.close_labware_latch()
            heater_shaker.set_and_wait_for_shake_speed(shaker)
            ctx.delay(minutes=0.5)
            heater_shaker.deactivate_shaker()    
            heater_shaker.open_labware_latch()
            ctx.move_labware(labware = reaction_plate,
                              new_location=mm_adapter,
                              use_gripper=True
                              )
            ctx.delay(minutes=isolation_time)

                # wash removal using 8 channel pipettes
            l = 1
            right_pipette.pick_up_tip()
            while l <= number_of_lane:
                right_pipette.flow_rate.aspirate = 10
                right_pipette.flow_rate.dispense = 200
                right_pipette.aspirate(transfer_amount, reaction_plate[f'A{l}'])
                right_pipette.dispense(transfer_amount, trash)
                #right_pipette.transfer(
                    #transfer_amount,
                    #reaction_plate[f'A{l}'],
                    #trash,
                    #new_tip='never'
                #)
                l += 1
            right_pipette.drop_tip()
        
        #left_pipette.drop_tip()
        #right_pipette.drop_tip()
        ctx.comment("This is the end of wash.")

    def deprotection(source_column = (1)):
        heater_shaker.open_labware_latch()
        #left_pipette.pick_up_tip()
        left_pipette.flow_rate.aspirate = 200
        left_pipette.flow_rate.dispense = 10
        
        for i in range(len(data[0])):
            if data[source_column][i] == 'SMAC_deprot':
                left_pipette.transfer(
                transfer_amount,
                tube_rack['A1'], #change later
                reaction_plate[f'{well_position[i]}'].bottom(3),
                blow_out=True,
                blowout_location='destination well'
                )
            elif data[source_column][i] == 'Fake_Wash':
                left_pipette.transfer(
                transfer_amount,
                tube_rack['A2'], #change later
                reaction_plate[f'{well_position[i]}'].bottom(3),
                blow_out=True,
                blowout_location='destination well'
                )
            else:
                print('Error')

        #left_pipette.drop_tip()
        ctx.move_labware(labware = reaction_plate,
                              new_location=hs_adapter,
                              use_gripper=True
                              )
        heater_shaker.close_labware_latch()
        heater_shaker.set_target_temperature(temperature)
        heater_shaker.set_and_wait_for_shake_speed(speed)
        ctx.delay(minutes=deprotection_time)  # 60 min
        heater_shaker.deactivate_heater()
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()

    def wash_deprotection(source_column: tuple = (1)):
        #removal of the solution with right pipette
        #source column[0] is number of lane
        
        #Removal of reaction mixture
        ctx.comment("This is the beginning of wash of PLP solution.")
        ctx.pause('Centrifuge the plate')
        right_pipette.flow_rate.aspirate = 10
        right_pipette.flow_rate.dispense = 200
        ctx.move_labware(labware = reaction_plate,
                              new_location=mm_adapter,
                              use_gripper=True
                              )
        ctx.delay(minutes=isolation_time)
        j = 1
        right_pipette.pick_up_tip()

        while j <= source_column:
            right_pipette.aspirate(transfer_amount, reaction_plate[f'A{j}'])
            right_pipette.dispense(transfer_amount, trash)
            j += 1
        right_pipette.drop_tip()

        #right_pipette.pick_up_tip()

        ### wash 2 times
        for i in range(2):
            k = 0
            right_pipette.pick_up_tip()
            while k <= number_of_lane - 1:
                right_pipette.flow_rate.aspirate = 200
                right_pipette.flow_rate.dispense = 10
                #left_pipette.pick_up_tip()
                right_pipette.transfer(
                    transfer_amount,
                    depro_wash_stock['A1'],
                    reaction_plate[f'A{k+1}'].bottom(3),
                    new_tip='never'
                )
            #left_pipette.mix(5, 200, magnetic_plate[f'{well_position[k]}'])
            #left_pipette.drop_tip()
                k += 1
            right_pipette.drop_tip()

            ctx.move_labware(labware = reaction_plate,
                                new_location=hs_adapter,
                                use_gripper=True
                                )
            heater_shaker.close_labware_latch()
            
            heater_shaker.set_and_wait_for_shake_speed(shaker)
            ctx.delay(minutes=0.5)
            heater_shaker.deactivate_shaker()    
            heater_shaker.open_labware_latch()
            ctx.move_labware(labware = reaction_plate,
                                new_location=mm_adapter,
                                use_gripper=True
                                )
            ctx.delay(minutes=isolation_time)

            # wash removal using 8 channel pipettes
            l = 1
            right_pipette.pick_up_tip()
            while l <= number_of_lane:
                right_pipette.flow_rate.aspirate = 10
                right_pipette.flow_rate.dispense = 200
                right_pipette.aspirate(transfer_amount, reaction_plate[f'A{l}'])
                right_pipette.dispense(transfer_amount, trash)
                l += 1
            right_pipette.drop_tip()

        ctx.comment("Equilibration of resin to buffer")
        #right_pipette.pick_up_tip()
        ### wash 2 times
        for i in range(2):
            right_pipette.pick_up_tip()
            k = 0
            while k <= number_of_lane - 1:
                right_pipette.flow_rate.aspirate = 200
                right_pipette.flow_rate.dispense = 10
                #left_pipette.pick_up_tip()
                right_pipette.transfer(
                    transfer_amount,
                    wash_stock['A1'],
                    reaction_plate[f'A{k+1}'].bottom(3),
                    new_tip='never'
                )
            #left_pipette.mix(5, 200, magnetic_plate[f'{well_position[k]}'])
            #left_pipette.drop_tip()
                k += 1
            right_pipette.drop_tip()

            ctx.move_labware(labware = reaction_plate,
                                new_location=hs_adapter,
                                use_gripper=True
                                )
            heater_shaker.close_labware_latch()
            
            heater_shaker.set_and_wait_for_shake_speed(shaker)
            ctx.delay(minutes=0.5)
            heater_shaker.deactivate_shaker()    
            heater_shaker.open_labware_latch()
            ctx.move_labware(labware = reaction_plate,
                                new_location=mm_adapter,
                                use_gripper=True
                                )
            ctx.delay(minutes=isolation_time)

            # wash removal using 8 channel pipettes
            l = 1
            right_pipette.pick_up_tip()
            while l <= number_of_lane:
                right_pipette.flow_rate.aspirate = 10
                right_pipette.flow_rate.dispense = 200
                right_pipette.aspirate(transfer_amount, reaction_plate[f'A{l}'])
                right_pipette.dispense(transfer_amount, trash)
                l += 1
            right_pipette.drop_tip()
        
        #equilibrating with 50 mM HEPES, 50 mM NaCl
        ctx.comment("This is the beginning of equilibrating method.")
        n = 0
        left_pipette.pick_up_tip()
        while n <= len(well_position) - 1:
            left_pipette.flow_rate.aspirate = 200
            left_pipette.flow_rate.dispense = 10
            left_pipette.transfer(
                transfer_amount,
                tube_rack['A2'],
                reaction_plate[f'{well_position[n]}'].bottom(3),
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            n += 1
        left_pipette.drop_tip()
        
        ctx.move_labware(labware = reaction_plate,
                                new_location=hs_adapter,
                                use_gripper=True
                                )
        heater_shaker.close_labware_latch()
            
        heater_shaker.set_and_wait_for_shake_speed(shaker)
        ctx.delay(minutes=0.5)
        heater_shaker.deactivate_shaker()    
        heater_shaker.open_labware_latch()
        ctx.move_labware(labware = reaction_plate,
                            new_location=mm_adapter,
                            use_gripper=True
                            )
        ctx.delay(minutes=isolation_time)

            # wash removal using 8 channel pipettes
        l = 1
        right_pipette.pick_up_tip()
        while l <= number_of_lane:
            right_pipette.flow_rate.aspirate = 10
            right_pipette.flow_rate.dispense = 200
            right_pipette.aspirate(transfer_amount, reaction_plate[f'A{l}'])
            right_pipette.dispense(transfer_amount, trash)
            l += 1
        right_pipette.drop_tip()


    def ending():
        ctx.comment("This is the beginning of ending method.")
        n = 0
        left_pipette.pick_up_tip()
        while n <= len(well_position) - 1:
            left_pipette.flow_rate.aspirate = 200
            left_pipette.flow_rate.dispense = 10
            left_pipette.transfer(
                transfer_amount,
                tube_rack['A2'],
                reaction_plate[f'{well_position[n]}'].bottom(3),
                blow_out=True,
                blowout_location='destination well',
                new_tip='never'
            )
            n += 1
        left_pipette.drop_tip()
        heater_shaker.open_labware_latch()
        ctx.comment("This is the end of ending method.")
    
    current_cycle = 1
    number_of_steps = ub_len -  acceptor
    
    ctx.pause('place the reaction plate to magnetic module')

    while current_cycle <= number_of_steps:
        #, PLP solution information
        deprotection(source_column=current_cycle-1) #, PLP solution information
        wash_deprotection(source_column=(number_of_lane))
        reaction(source_column=(number_of_lane, current_cycle))
        wash_reaction(source_column=(number_of_lane, current_cycle))
        current_cycle += 1
    
    ending()
    