-- Simulation Settings
-- Note all units are metric

-- System state
l1 = 0.5
l2 = 0.5
m1 = 0.2
m2 = 0.2
pixelsPerMeter = 300

-- Initial Conditions
dTheta1 = 0
dTheta2 = 0
Theta1 = 3.14159/2
Theta2 = 3.14159/2

g = 9.80665

file = 0

-- f1 and f2 are the second derivatives of theta1 and theta2. They'll be called a lot and so are defined here for later use
function f1(theta1, theta2, dtheta1, dtheta2)
	numerator = g*m1*math.sin(theta1)+0.5*g*m2*math.sin(theta1-2*theta2)+0.5*g*m2*math.sin(theta1)+0.5*l1*m2*math.sin(2*(theta1-theta2))*math.pow(dtheta1, 2)+l2*m2*math.sin(theta1-theta2)*math.pow(dtheta2, 2)
	denominator = -l1*(m1-m2*math.pow(math.cos(theta1-theta2), 2)+m2)
	return numerator/denominator
end
function f2(theta1, theta2, dtheta1, dtheta2)
	mu = 1+m1/m2
	numerator = g*mu*(math.sin(theta1)*math.cos(theta1-theta2)-math.sin(theta2))+(mu*l1*math.pow(dtheta1, 2) + l2*math.pow(dtheta2, 2)*math.cos(theta1-theta2))*math.sin(theta1-theta2)
	denominator = l2*(mu-math.pow(math.cos(theta1-theta2), 2))
	return numerator/denominator
end

function love.load()
	love.graphics.setBackgroundColor(104, 136, 248) --set the background color to a nice blue
    love.window.setMode(900, 700) --set the window dimensions to 650 by 650

    -- I output the state of the system to a text file so the motion can be graphed and analysed elsewhere
	file = love.filesystem.newFile("data.txt")
    file:open("w")
    file:write("theta1, theta2, dtheta1, dtheta2, pot1, pot2, kin1, kin2 \r\n")
end

function love.quit()
	file:close()
end

-- This is where the meat of the solving happens. Note this is a 4th order Runga Kutta solver
function love.update(dt)

	-- The slope is calculated at various points along the timestep of the graph. And note, we have to find these slopes for acceleration and then velocity
	-- to get back to displacement.

	L1_1 = dTheta1*dt -- Slope of angular displacement (angular velocity), used for calculating the next angular displacement
	L1_2 = dTheta2*dt
	K1_1 = f1(Theta1, Theta2, dTheta1, dTheta2)*dt -- Slope of angular velocity (angular acceleration), used for calculating next angular velocity
	K1_2 = f2(Theta1, Theta2, dTheta1, dTheta2)*dt

	L2_1 = (dTheta1+K1_1/2)*dt
	L2_2 = (dTheta2+K1_2/2)*dt
	K2_1 = f1(Theta1+L1_1/2, Theta2+L1_2/2, dTheta1+K1_1/2, dTheta2+K1_2/2)*dt
	K2_2 = f2(Theta1+L1_1/2, Theta2+L1_2/2, dTheta1+K1_1/2, dTheta2+K1_2/2)*dt


	L3_1 = (dTheta1+K2_1/2)*dt
	L3_2 = (dTheta2+K2_2/2)*dt
	K3_1 = f1(Theta1+L2_1/2, Theta2+L2_2/2, dTheta1+K2_1/2, dTheta2+K2_2/2)*dt
	K3_2 = f2(Theta1+L2_1/2, Theta2+L2_2/2, dTheta1+K2_1/2, dTheta2+K2_2/2)*dt

	L4_1 = (dTheta1+K3_1)*dt
	L4_2 = (dTheta2+K3_2)*dt
	K4_1 = f1(Theta1+L3_1, Theta2+L3_2, dTheta1+K3_1, dTheta2+K3_1)*dt
	K4_2 = f2(Theta1+L3_1, Theta2+L3_2, dTheta1+K3_1, dTheta2+K3_2)*dt


	-- Once we have all of our slopes we average the slopes to find the best slope over the whole timestep. This is added
	-- to velocity and then displacement to get our final solution for this timestep
	dTheta1 = dTheta1 + (K1_1+2*K2_1+2*K3_1+K4_1)/6
	dTheta2 = dTheta2 + (K1_2+2*K2_2+2*K3_2+K4_2)/6
	
	Theta1 = Theta1 + (L1_1+2*L2_1+2*L3_1+L4_1)/6
	Theta2 = Theta2 + (L1_2+2*L2_2+2*L3_2+L4_2)/6

	-- Here I calculate the potential energy and kinetic energy expressly for the purpose of outputing the data and graphing it later
	pot1 = -g*l1*m1*math.cos(Theta1)
	pot2 = g*m2*(-l1*math.cos(Theta1) - l2*math.cos(Theta2))
	kin1 = math.pow(l1, 2)*m1*math.pow(dTheta1, 2)/2
	kin2 = m2*(math.pow(l1,2)*math.pow(dTheta1,2) + 2*l1*l2*math.cos(Theta1 - Theta2)*dTheta1*dTheta2 + math.pow(l2, 2)*math.pow(dTheta2,2))/2

	file:write(Theta1..", "..Theta2..", "..dTheta1..", "..dTheta2..", "..pot1..", "..pot2..", "..kin1..", "..kin2.." \r\n")
end

function love.draw()
	-- System Anchor
	x0 = 450
	y0 = 350
	-- We have  to convert from our coordinate system back into cartesian coordinates for Love2D
	-- Position of Bob 1
	x1 = x0 + l1*math.sin(Theta1) * pixelsPerMeter
	y1 = y0 + l1*math.cos(Theta1) * pixelsPerMeter
	-- Position of Bob 2
	x2 = x1 + l2*math.sin(Theta2) * pixelsPerMeter
	y2 = y1 + l2*math.cos(Theta2) * pixelsPerMeter
	
	-- Draw the pendulum strings
	love.graphics.setColor(0, 0, 0)
	love.graphics.line(x0, y0, x1, y1)
	love.graphics.line(x1, y1, x2, y2)

	-- Draw the pendulum bobs
	love.graphics.setColor(193, 47, 14)
    love.graphics.circle("fill", x1, y1, 10) -- First Bob 1
    love.graphics.circle("fill", x2, y2, 10) -- Second Bob 2
end