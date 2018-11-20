# Jiuyang Bai 1D following v1.1 Carrot3
# Data collected for first project
# Carrot3 experiment studying speed control in one dimensional follwing

# The only difference from Carrot2 are the conditions: d0 = 1/4/8  v0 = 0.8/1.2  dv = -0.3/0/0.3  a = 1

# 08/09/2017

import math
import lights # turns on the lights
import random
import csv
#from os.path import exists
import os

# Vizard Imports
import viz
import vizact
import viztracker


# The following libraries are won't work outside of the VENLab without their 
# associated dependencies, but are requried for experiments within the VENLab.
# When experiment is ready to be moved to the VENLab, we'll re-include them.

import oculus
import emergencyWalls

#####################################################################################
# Constants


IPD = viz.input('Please enter IPD:')

# Set to True when ready to have the experiment write data
DATA_COLLECT = True
 
# If program should run practice trials
DO_PRACTICE = False

# If run free walk trials
DO_FREEWALK = False

# If program crashes, start trials here
START_ON_TRIAL = 1

# Total number of trials in experiment
FREEWALK_TRIALS = 4 # 4 trials each session
FREEWALK_SESSIONS = 2 # 1 session before practice 1 session after experiment
PRACTICE_TRIALS = 4
TOTAL_TRIALS = 90    # 3(d0) * 2(v0) *  3(dv) * 5(reps) = 90
TRIAL_LENGTH = 12 # 12 seconds


# Used for file naming, currently placeholders
EXPERIMENT = '1Dfollowing_v1.1_Carrot3'
EXPERIMENTER_NAME = 'jiuyangbai'

# inputFile = 'Data/jiuyangbai/1Dfollowing_v1.1_Carrot2/Subject_00/Input/trial001.csv'


# Orientation constants
POLE_TRIGGER_RADIUS = 0.3 # How close participant must be to home pole
THRESHOLD_THETA = 15 # Maximum angle participant can deviate when looking at orienting pole
ORIENT_TIME = 3 # How long participant must orient onto pole

# The dimension of the room space used for experiment
DIMENSION_X = 9.0 # the length of the shorter side in meter
DIMENSION_Z = 11.0 # the lengtth of the longer side in meter
DIAGONAL = (DIMENSION_X**2 + DIMENSION_Z**2)**(1.0/2)# The length of the diagonal line of the experimental space

ROOM_ANGLE = math.atan(DIMENSION_X/DIMENSION_Z) # the anger (in radian) between the diagonal and the shorter edge of the room

# Home and Orient Pole positions (x,z,y)
HOME_POLE = [[DIMENSION_X/2, 0.0, DIMENSION_Z/2], [-DIMENSION_X/2, 0.0, -DIMENSION_Z/2]]
ORI_POLE = [[DIMENSION_X/2 - DIMENSION_X/3, 0.0, DIMENSION_Z/2 - DIMENSION_Z/3], \
			[-DIMENSION_X/2 + DIMENSION_X/3, 0.0, -DIMENSION_Z/2 + DIMENSION_Z/3]]


# Describe the end-trial trigger line (end line) in the intersense coordinate system
# the end line is perpendicular to walking direction
K = -DIMENSION_X/DIMENSION_Z # The slope of the end line
END_DIS = 2.0 # The distance between end line to home pole position of following trial
# The intercept of the end line for two home poles respectively
B = [-(DIAGONAL/2 - END_DIS) / math.cos(ROOM_ANGLE), (DIAGONAL/2 - END_DIS) / math.cos(ROOM_ANGLE)]

#####################################################################################
# Control Options
# Sets the controls / displays
# Oculus Rift Selection unlikely to work outside the VENLab


# Dialog box asking for type of control
OCULUS = 'Oculus Rift'
MONITOR = 'PC Monitor'
controlOptions = [OCULUS,MONITOR]
controlType = controlOptions[viz.choose('How would you like to explore? ', controlOptions)]



# Use keyboard controls
# Controls:
# q - Strafe L		w - Forward		e - Strafe R
# a - Turn L		s - Back		d - Turn R
#
# y - Face Up		r - Fly Up
# h - Face Down		f - Fly Down
if controlType == MONITOR:
	headTrack = viztracker.Keyboard6DOF()
	link = viz.link(headTrack, viz.MainView)
	headTrack.eyeheight(1.6)
	link.setEnabled(True)
	viz.go()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
# Use Oculus Rift
elif controlType == OCULUS:
	# viz.fullscreen.x is the inital x position when run
	viz.ipd(IPD)	

	#viz.vsync(0)
	#viz.setDisplayMode(0,0,0,60)
	viz.go()
	
	# add intersense tracker
	isense = viz.add('intersense.dle')
	ISTracker = isense.addTracker(port=5001,station=3)
	ISTracker.setEnhancement(2) 
	ISTracker.setSensitivity(4)
	ISTracker.setShockSuppression(2)
	ISTracker.setAccelSensitivity(4)	
	# add Oculus tracker
	OVRTracker = oculus.Rift().getSensor()
	
	# add the virtual tracker, link it to the MainView and set an offset to get the eye position
	virtual = viz.addGroup()
	link = viz.link(virtual, viz.MainView)
	link.preTrans([0, -0.055, -0.073])
	
	# tracker initialization flag
	init = False
	
	# yaw correction flag
	corrected = False

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Use Oculus Rift
#elif controlType == OCULUS:
#	hmd = oculus.Rift()
#	viz.setOption('viz.fullscreen.x','0')
#	viz.setOption('viz.fullscreen.y','0')
#	vizconnect.go('vrpnOculusMerge.pyc')
#	hmdName = 'oculus'

viz.clip(.001,1000) # Adjust size of the viewing frustum



######################################################################################################
# Experiment Loop Helper Functions
def goToStage(nextTrialStage):
	global trial_stage
	trial_stage = nextTrialStage
	print 'Going to: ' + trial_stage

	
def endLine(x):
	# takes the x value and gives the z value of the end trigger line defined in the intersense coordinate system
	return x*K + B[trial_num%2]
	
def moveTarget(traj, index):
	
	global posIndex, time, targetTraj
	
	if trial_num%2 == 1: 
		# start from (-, 0, -) corner of the room	
		dz = float(traj[index][1])*math.cos(ROOM_ANGLE) - float(traj[index][0])*math.sin(ROOM_ANGLE)
		dx = float(traj[index][1])*math.sin(ROOM_ANGLE) + float(traj[index][0])*math.cos(ROOM_ANGLE)
	else:                
		# start from (+, 0, +) corner of the room
		dz = -float(traj[index][1])*math.cos(ROOM_ANGLE) + float(traj[index][0])*math.sin(ROOM_ANGLE)
		dx = -float(traj[index][1])*math.sin(ROOM_ANGLE) - float(traj[index][0])*math.cos(ROOM_ANGLE)
	
	models['targetPole'].setPosition([HOME_POLE[trial_num%2][0]+dx, 0, HOME_POLE[trial_num%2][2]+dz])
	

def countDown(t):
	global time, stamped, time_stamp
	timesUp = False
	if not stamped:
		time_stamp = time
		stamped = True
	if time - time_stamp > t:
		timesUp = True
	return timesUp
	
def writeCSVFile(fileName, data, time):
	strData = [str(round(t,4)) for t in data+[time]]
	file = open(fileName, 'a')
	file.write(','.join(strData)+'\n')
	file.close()
	
def relativeOrientation(pos1, pos2):
	xrel = round(pos2[0]-pos1[0],4)
	zrel = round(pos2[2]-pos1[2],4)
	theta = 0
	if zrel == 0.0 and xrel > 0:
		theta = math.pi/2
	elif zrel == 0.0:
		theta = math.pi/2*3
	else:
		theta = math.atan(round(xrel,4)/round(zrel,4))
		if zrel < 0:
			theta += math.pi
		if zrel > 0 and xrel < 0:
			theta += math.pi*2
	return theta
	
def facing(lookingObjectPosn, lookedAtObjectPosn, lookingObjectYaw, thresholdTheta):
	"""lookingObjectPosn: position of object that is looking
	lookedAtObjectPosn: position of object looked at
	lookingObjectYaw: yaw of the object that is looking (degrees)
	thresholdTheta: viewing angle must be +/- this amount in order to be considered 'looking at' the object. degrees

	return: bool, whether the looking object is facing the looked-at object
	>>> universals.facing([0,0,0],[1,0,5],0,20)
	True
	>>> universals.facing([3,0,3],[1,0,0],210,20)
	True
	"""
	degRelOrientation = 180.0/math.pi*relativeOrientation(lookingObjectPosn, lookedAtObjectPosn) #radians
	degRelOrientation = (degRelOrientation+180)%360-180
	return math.fabs(degRelOrientation-lookingObjectYaw)<thresholdTheta
	
def distance(x,y,a,b):
	return ((x-a)**2+(y-b)**2)**.5
	
def inRadius(pos, center, radius):
	#This method takes in two poadsitions in [x,y,z] form and then returns true if the distance between them is less than the radius given
	if pos == '' or center == '' or radius == '':
		return False
	return (distance(pos[0],pos[2],center[0],center[2]) <= radius)
	

	
#######################################################################################################
# Experiment

subject = viz.input('Please enter the subject number:','')
if subject < 10:
	subject = '0' + str(subject)

# Set output directory for writing data
# Location: Same folder as this .py file/Data/dachner/[experiment_name]
output_dir = '/'.join(['Data',EXPERIMENTER_NAME,EXPERIMENT,'Subject_' + str(subject),'Output'])
input_dir = '/'.join(['Data',EXPERIMENTER_NAME,EXPERIMENT,'Subject_' + str(subject),'Input'])


# loads condition data
inputFile = input_dir + '/conditions.csv'
with open(inputFile, 'rb') as data:
	rows = csv.reader(data)
	conditions = [[value for key, value in enumerate(row)] for row in rows]
condition = ''

viz.clearcolor(0,0.4,1.0) # blue world
models = {}
models['homePole'] = viz.add('Models/bluePole.3ds')
models['orientPole'] = viz.add('Models/redPole.3ds')
models['targetPole'] = viz.add('Models/greenPole.3ds')
models['ground'] = viz.add('Models/ground4.osgb')
# Adjust models size
models['homePole'].setScale([0.6,0.45,0.6]) # the original size = [0.4m 3m 0.4m]
models['targetPole'].setScale([1.0,0.66,1.0]) # the original size = [0.4m 3m 0.4m]
# the transparency 
_alpha = 0.0
# Hide loaded models
models['homePole'].visible(viz.OFF)
models['orientPole'].visible(viz.OFF)
models['targetPole'].visible(viz.OFF)


# Sounds
sounds={}
sounds['intro'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/intro.mp3')
sounds['freeTarget'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/freeTarget.mp3')
sounds['practiceHome'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/practiceHome.mp3')
sounds['testHome'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/testHome.mp3')
sounds['postFree'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/postFree.mp3')
sounds['begin'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/begin.mp3')
sounds['end'] = viz.addAudio('Sounds/jiuyangbai/1Dfollowing_v1.1_Carrot3/end.mp3')


# Initial data_collection, regardless of DATA_COLLECT
data_collect = False

# Initializes trial_stage, which is the current step of a given trial
# First part is the name of the specific stage (pretrial)
# Second part is practice (01) or experimental trials (02)
# Third part is the stage's position in the order of all stages (01)
goToStage('pretrial_00_01')

# Initializa setting to start free walk trials
if DO_FREEWALK == True:
	is_freewalk = True
else:
	is_freewalk = False
	goToStage('pretrial_01_01')

# Initially setting to start practice trials
if DO_PRACTICE == True:
	is_practice = True
else:
	is_practice = False
	goToStage('pretrial_02_01')
	
# Starts trial count at 1, unless otherwise specifed
if START_ON_TRIAL > 1:
	trial_num = START_ON_TRIAL
	
else:
	trial_num = 1
freewalk_session = 1
# Time counter set to 0
time = 0

instruction = True
stamped = False

screenshot = 1
videorecord = 1

########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################

def masterLoop(num):
	# global variables within masterLoop
	global DATA_COLLECT, DO_PRACTICE, is_practice, is_freewalk, data_collect, trial_stage, trial_num, \
	freewalk_session, time, time_stamp, cur_pos, posIndex, conditions, condition, K, B, targetTraj, _alpha, \
	stamped, controlType,instruction, screenshot, videorecord, \
	old_OVROri, old_camOri, init, corrected, virtual #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# Time elapsed since the last run of masterLoop and then added to the global time
	frame_elapsed = viz.getFrameElapsed()
	time += frame_elapsed
	
	viz.setOption('viz.AVIRecorder.maxWidth','1024')
	viz.setOption('viz.AVIRecorder.maxHeight','720')
	viz.setOption('viz.AVIRecorder.fps','30')
	

	if os.path.isfile(output_dir + '/image'+ str(screenshot) +'.bmp') == True:
		screenshot += 1
	if os.path.isfile(output_dir + '/video'+ str(videorecord) +'.avi') == True:
		videorecord += 1
		
	vizact.onkeydown('p', viz.window.screenCapture, output_dir + '/image'+ str(screenshot) +'.bmp')
	vizact.onkeydown('b', viz.window.startRecording, output_dir + '/video'+ str(videorecord) +'.avi')
	vizact.onkeydown('e', viz.window.stopRecording)
	

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	if controlType == OCULUS:
		if not init:
			try:
				old_OVROri = OVRTracker.getEuler()
				# initial the camera euler by Interesense (Actually just want yaw. Pitch and roll will be replaced by Oculus later)
				old_camOri = ISTracker.getEuler()
				# replace pitch and roll by Oculus because its more sensitive and doesn't give too much error
				old_camOri[1] = old_OVROri[1]
				old_camOri[2] = old_OVROri[2]
			except NameError:
				print('tracker haven''t been initialized')
			else:
				init = True

		else:
		
			# get the current orientation from trackers
			ISOri = ISTracker.getEuler()
			OVROri = OVRTracker.getEuler()
			
			# apply the change of Oculus orientation to the camera
			d_OVROri = [OVROri[0] - old_OVROri[0],\
						OVROri[1] - old_OVROri[1],\
						OVROri[2] - old_OVROri[2]]
						
			camOri = [old_camOri[0] + d_OVROri[0],\
					  old_camOri[1] + d_OVROri[1],\
					  old_camOri[2] + d_OVROri[2]]


			# Update MainView by intersense position and intersense-corrected Oculus orientation
			virtual.setPosition(ISTracker.getPosition())
			virtual.setEuler(camOri)
			# The orientation at this frame become the old orientation for next frame
			old_camOri = camOri
			old_OVROri = OVROri

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	# Current position and roation of the participant
	cur_pos = viz.get(viz.HEAD_POS)
	cur_rot = viz.get(viz.HEAD_ORI)
	
	#>> Will only work in VENLab <<
	emergencyWalls.popWalls(cur_pos) # Pops up the Emergency Walls when participant is close to physical room edge.
	
	
	
	
	
	
#	some testing code
#	data = [1,1,1,1]
#	# the format of condition is "d0_v0_dv"
#	writeCSVFile(output_dir + '/framerate_test.txt', data, viz.getFrameElapsed())
#	print('haha')
	
	
	
	
	
	##################
	# Begin Freewalk #
	##################
	
	if is_freewalk == True:
		
		
		# Writes Position and Rotation, but only when DATA_COLLECT set to True
		# and trial_stage has set data_collect to True
		if DATA_COLLECT and data_collect:

			
			target_loc = models['targetPole'].getPosition()

			# File Name: subj_#_condition_#_freewalk_#_position.csv
			# Format: Target_x; Target_y; Participant_x; Participant_y; trial time
			
			data = [target_loc[0],target_loc[2],cur_pos[0], cur_pos[2]]
			writeCSVFile(output_dir + '/' + 'position_' + 'freewalk_' + 's' + str(freewalk_session) + \
			'_subj_' + str(subject) + '_trial_' + str(trial_num) + '.csv', data, time)
			
			# Orientation
			data = [cur_rot[0], cur_rot[1], cur_rot[2]]
			writeCSVFile(output_dir + '/' + 'orientation_' + 'freewalk_' + 's' + str(freewalk_session) + \
			'_subj_' + str(subject) + '_trial_' + str(trial_num) + '.csv', data, time)
		
		#########
		# 00 01 Freewalk Pretrial: sets up practice trial, establishes pole locations
		if trial_stage == 'pretrial_00_01':

			print '> Start Free Walking Session ' + str(freewalk_session) + ' Trial ' + str(trial_num) + ' ----------------------------------'
		
			# Set position of home pole (where participant stands to start trial)			
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				models['homePole'].alpha(1.0)
				models['homePole'].visible(viz.ON)

			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
			
			if trial_num == 1:
				if freewalk_session == 1:
					sounds['intro'].play()
				elif freewalk_session == 2:
					sounds['postFree'].play()
				
			
			# Move to Stage 2
			goToStage('orient_00_02')
		
			
		#########
		# 00 02 Orienting to Pole: Give time for participant to orient to the pole
		elif (trial_stage == 'orient_00_02'):
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if controlType == OCULUS:
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					
					# Yaw correction
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					
			if (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
			and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				
				# Move to stage 3
				goToStage('wait for orientation')
				corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
				# Current time
				time_stamp = time
		
		#########
		# wait for orientation
		elif (trial_stage == 'wait for orientation'):
			if countDown(ORIENT_TIME):
				goToStage('inposition_00_03')
			if not (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS) \
			and facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				stamped = False
		
		#########
		# 00 03 Freewalk In Position: proceeds once participant is standing on home and facing orient for three seconds
		elif (trial_stage == 'inposition_00_03') and (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
		and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)) \
		and (time-time_stamp >= ORIENT_TIME):
			
			print 'Free walk start'
			
			# Turn off home pole and orient pole
			models['homePole'].visible(viz.OFF)			
			models['orientPole'].visible(viz.OFF)
			models['targetPole'].setPosition(HOME_POLE[(trial_num+1)%2])
			models['targetPole'].visible(viz.ON)
			if trial_num > 1:
				sounds['begin'].play()
			
			if trial_num == 1 and freewalk_session == 1:				
				sounds['freeTarget'].play()
				
			# Start to collect data
			data_collect = True
		
			# Move to Stage 4
			goToStage('target_00_04')
			time = 0
			
		#########
		# 00 04 Freewalk: Participants Moves
		elif (trial_stage == 'target_00_04'):
		
				
			# Detects participant location, moves to Stage 5 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])):
					
				goToStage('endtrial_00_05')
		
		#########
		# 00 05 End Freewalk Trial: Close out the trial and reset values for next practice trial or start Experiment
		elif trial_stage == 'endtrial_00_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			print 'End Freewalk Trial ' + str(trial_num)
			data_collect = False


			# End Check: When trial_num is greater than FREEWALK_TRIALS, end practice and start experiment block
			if trial_num == FREEWALK_TRIALS:
				print '>> End Freewalk Session<<'
	
				if freewalk_session >= 2:
					print '>>> End Experiment <<<'
					goToStage('NULL')
					if instruction:						
						sounds['end'].play()
						instruction = False
				else:
					goToStage('pretrial_01_01')
					is_freewalk = False
					trial_num = 1
					freewalk_session += 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_00_01')				
				instruction = True
		
			
			
	##################
	# Begin practice #
	##################
	
	elif is_practice == True:
		
	
		#########
		# 01 01 Practice Pretrial: sets up practice trial, establishes pole locations
		if trial_stage == 'pretrial_01_01':

			print '> Start Practice Trial ' + str(trial_num) + ' ----------------------------------'
			
			# Loads input files

			inputFile = input_dir + '/practice' + str(trial_num) + '.csv'

			with open(inputFile, 'rb') as data:
				rows = csv.reader(data)
				targetTraj = [[value for key, value in enumerate(row)] for row in rows]
	
			# will be used later as the index of targetTraj to update target pole position
			posIndex = 0
	
			# Move to Stage 2
			goToStage('wait for instruction')
		
			
			
		#########
		# wait for instruction
		elif (trial_stage == 'wait for instruction'):
			if trial_num == 1: 
				if instruction:
					sounds['practiceHome'].play()
					instruction = False

				if countDown(1):
					goToStage('orient_01_02')
					stamped = False
				
			else:
				goToStage('orient_01_02')
			
			
			
		#########
		# 01 02 Orienting to Pole: Give time for participant to orient to the pole
		elif (trial_stage == 'orient_01_02'):

			# Set position of home pole (where participant stands to start trial)
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				
			if _alpha < 1.0:
				models['homePole'].alpha(_alpha)
				models['homePole'].visible(viz.ON)
				_alpha += 0.016
	
			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
		
			
				
				
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					# Yaw correction			
			if controlType == OCULUS:	
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					
			if (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
			and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				
				# Move to stage 3
				goToStage('wait for orientation')
				corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
				_alpha = 0.0
				# Current time
				time_stamp = time
				
		#########
		# wait for orientation
		elif (trial_stage == 'wait for orientation'):
			if countDown(ORIENT_TIME):
				goToStage('inposition_01_03')
			if not (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS) \
			and facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				stamped = False
				
				
		#########
		# 01 03 Practice In Position: proceeds once participant is standing on home and facing orient for three seconds
		elif (trial_stage == 'inposition_01_03') and (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
		and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)) \
		and (time-time_stamp >= ORIENT_TIME):
			
			print 'Practice Target Appears'
			
			# Turn off home pole and orient pole
			models['homePole'].visible(viz.OFF)
			models['orientPole'].visible(viz.OFF)


			# Initializes target pole at home location, which is the reference point for moveTarget function to update target position)
			models['targetPole'].setPosition(HOME_POLE[trial_num%2]) 
			moveTarget(targetTraj, 0) # update the position once to create initial distance before it appears
			models['targetPole'].visible(viz.ON) # target pole appears
			sounds['begin'].play()
			
			# Move to Stage 4
			goToStage('target_01_04')
			time = 0
			
	
		
		#########
		# 01 04 Moving Pole: Target Moves
		elif (trial_stage == 'target_01_04'):

			# Target Moves
			moveTarget(targetTraj, posIndex)
			posIndex = int(time/TRIAL_LENGTH * len(targetTraj))
			
				
			# Detects participant location, moves to Stage 6 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])) or \
				posIndex >= len(targetTraj):
				goToStage('endtrial_01_05')
		
		#########
		# 01 05 End Practice Trial: Close out the trial and reset values for next practice trial or start Experiment
		elif trial_stage == 'endtrial_01_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			print 'End Practice Trial ' + str(trial_num)
			

			# End Check: When trial_num is greater than PRACTICE_TRIALS, end practice and start experiment block
			if trial_num >= PRACTICE_TRIALS:
				print '>> End Practice <<'
				goToStage('pretrial_02_01')
				is_practice = False
				trial_num = 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_01_01')
				instruction = True
		
			time = 0
			
			
	####################
	# Begin Experiment #
	####################
	
	elif is_practice == False:
		
		
		condition = ', '.join([conditions[trial_num][1], conditions[trial_num][2], conditions[trial_num][3]])
		
		# Writes Position and Rotation, but only when DATA_COLLECT set to True
		# and trial_stage has set data_collect to True
		if DATA_COLLECT and data_collect:
			
			
			# Location of Target Pole
			target_loc = models['targetPole'].getPosition()
			
			# Format: Target_x; Target_y; Participant_x; Participant_y; trial time

			
			data = [target_loc[0],target_loc[2],cur_pos[0], cur_pos[2]]
#			data = [0.0987,1.1235,2.1345,0.1234]
			# the format of condition is "d0_v0_dv"			
			writeCSVFile(output_dir + '/' + 'position_' + condition + '_subj_' + str(subject) + '_trial_' + str(trial_num) + \
			'.csv', data, time)
			
			# Orientation
			data = [cur_rot[0], cur_rot[1], cur_rot[2]]
#			data = [0.0987,0.1232,1.1645,0.1234]
			writeCSVFile(output_dir + '/' + 'orientation_' + condition + '_subj_' + str(subject) + '_trial_' + str(trial_num) + \
			'.csv', data, time)
			
			# log IPD
			if trial_num == 1:
				file = open(output_dir + '/' + 'subj_' + str(subject) + \
				'_IPD_' + str(IPD) + '.txt', 'a')
				file.write('IPD = ' + str(IPD))
				file.close()
			
			
		#########
		# 02 01 Experiment Pretrial: sets up trial, establishes pole locations
		if trial_stage == 'pretrial_02_01':
			
			
			# Print start of trial, trial #, and type of trial [pos][speed][turn]
			print '> Start Trial ' + str(trial_num) + ': ' + condition + ' ----------------------------------'
					
			# Loads input files

			if trial_num < 10:
				inputFile = input_dir + '/trial00' + str(trial_num) + '.csv'
			else:
				inputFile = input_dir + '/trial0' + str(trial_num) + '.csv'
					
			with open(inputFile, 'rb') as data:
				rows = csv.reader(data)
				targetTraj = [[value for key, value in enumerate(row)] for row in rows]
	
			# will be used later as the index of targetTraj to update target pole position
			posIndex = 0
			
			
			# Move to Stage 2
			goToStage('wait for instruction')
		
		#########
		# wait for instruction
		elif (trial_stage == 'wait for instruction'):
			if trial_num == 1: 
				if instruction:
					sounds['testHome'].play()
					instruction = False

				if countDown(1):
					goToStage('orient_02_02')
					stamped = False
				
			else:
				goToStage('orient_02_02')
				
				
		#########
		# 02 02 Orienting to Pole: Give time for participant to orient to the pole
		elif (trial_stage == 'orient_02_02'):
			
			## Set position of home pole (where participant stands to start trial)
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				
			if _alpha < 1.0:
				models['homePole'].alpha(_alpha)
				models['homePole'].visible(viz.ON)
				_alpha += 0.016
	
			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
			
			
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
			if controlType == OCULUS:
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					# Yaw correction
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			if (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
			and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				
				
				# Move to stage 3
				goToStage('wait for orientation')
				corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
				_alpha = 0.0
				# Current time
				time_stamp = time
		
			
		#########
		# wait for orientation
		elif (trial_stage == 'wait for orientation'):
			if countDown(ORIENT_TIME):
				goToStage('inposition_02_03')
			if not (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS) \
			and facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)):
				stamped = False	
	
		
		#########
		# 02 03 In Position: proceeds once participant is standing on home and facing orient
		elif (trial_stage == 'inposition_02_03') and (inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS)) \
		and (facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)) \
		and (time-time_stamp >= ORIENT_TIME):

			# Turn off home and orient poles
			models['homePole'].visible(viz.OFF)
			models['orientPole'].visible(viz.OFF)
				
			# Initializes target pole at home location, which is the reference point for moveTarget function to update target position)
			models['targetPole'].setPosition(HOME_POLE[trial_num%2]) 
			moveTarget(targetTraj, 0) # update the position once to create initial distance before it appears
			models['targetPole'].visible(viz.ON) # target pole appears
			sounds['begin'].play()
			
			# Turn on data collection for this trial
			data_collect = True
				
			print 'Target Appears'
				
			# Move to Stage 5
			goToStage('target_02_04')
			
			time = 0

			
		#########
		# 02 04 Moving Pole: Target Moves
		elif (trial_stage == 'target_02_04'):
			
			# Target Moves			
			moveTarget(targetTraj, posIndex)
			posIndex = int(time/TRIAL_LENGTH * len(targetTraj))
			
			# Detects participant location, moves to Stage 6 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])) or \
				posIndex >= len(targetTraj):
				goToStage('endtrial_02_05')
				
			#########
		# 02 05 End Trial: Close out the trial and reset values for next trial
		elif trial_stage == 'endtrial_02_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			# End data collection for this trial
			data_collect = False
			
			print 'End Trial ' + str(trial_num)
			
						
			
			# When trial_num is greater than TOTAL_TRIALS, end experiment
			if trial_num >= TOTAL_TRIALS:
				is_freewalk = True
				goToStage('pretrial_00_01')
				trial_num = 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_02_01')
				instruction = True

	

# Restarts the loop, at a rate of 60Hz
viz.callback(viz.TIMER_EVENT,masterLoop)
viz.starttimer(0,1.0/90,viz.FOREVER)

