# Python script for numerically modeling the expansion of the universe.
# Written by Nathan Reed, August 2011; released into the public domain.

import math

# Constants

hubbleConst = 7.2e-11				# in inverse years
gravConst = 66450					# in m^3 / kg / yr^2
matterDensity = 2.53e-27			# in kg / m^3
radiationDensity = 5.6e-31			#  ...
darkEnergyDensity = 6.78e-27		#  ...

darkEnergyW = -1.0					# dimensionless equation-of-state parameter
darkEnergyFactor = 3.0 * darkEnergyW + 1.0
darkEnergyPower = 3.0 * darkEnergyW + 2.0

timeStep = 1.0e7					# in years

# Run simulation backward from present day until scale factor is < 0.01

lResultsBackward = []

curTime = 0.0
curScale = 1.0
curHubble = hubbleConst

while curScale > 0.01:
	matterTerm = matterDensity / (curScale ** 2)
	radiationTerm = 2.0 * radiationDensity / (curScale ** 3)
	darkEnergyTerm = darkEnergyFactor * darkEnergyDensity / (curScale ** darkEnergyPower)
	
	accel = -4.0/3.0 * math.pi * gravConst * (matterTerm + radiationTerm + darkEnergyTerm)
	curHubble = curHubble - timeStep * accel
	curScale = curScale - timeStep * curHubble
	curTime = curTime - timeStep

	lResultsBackward.append((curTime / 1.0e9, curScale))

# Run simulation forward from present day for another 10 billion years

lResultsForward = []

curTime = 0.0
curScale = 1.0
curHubble = hubbleConst

while curTime < 1.0e10:
	matterTerm = matterDensity / (curScale ** 2)
	radiationTerm = 2.0 * radiationDensity / (curScale ** 3)
	darkEnergyTerm = darkEnergyFactor * darkEnergyDensity / (curScale ** darkEnergyPower)
	
	accel = -4.0/3.0 * math.pi * gravConst * (matterTerm + radiationTerm + darkEnergyTerm)
	curHubble = curHubble + timeStep * accel
	curScale = curScale + timeStep * curHubble
	curTime = curTime + timeStep

	lResultsForward.append((curTime / 1.0e9, curScale))

# Write results to a tab-delimited text file

lResultsBackward.reverse()
lResults = lResultsBackward + [(0.0, 1.0)] + lResultsForward
f = open('bigbang.txt', 'w')
for result in lResults:
	f.write('%f\t%f\n' % result)

print('Results written to bigbang.txt.')
