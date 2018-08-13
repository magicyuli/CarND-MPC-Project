[update_equation]: ./img/update_eq.png "Update equation"

# The Model

### state

The state includes vehicle coordinates `x` and `y`, velocity `v`, vehicle angle `psi`, cross track error `cte`
and error in vehicle angle relative to the reference trajectory in the vehicle's coordinate.

### actuators
Actuators include acceleration `a` and the angular acceleration `delta`.

### update equations
The update equation is as follows:
![Update equation][update_equation]

# Timestep Length and Elapsed Duration (N & dt)
My final decision on `N` and `dt` are `30` and `0.05`. The reasoning behind it is when the car travels at a velocity of 
`20 mph`, it travels almost `1 m` per `0.1 s`, and we should at least check on the controls every `1 m`. `dt=0.05` provides 
finer control. `N=30` allows to see far enough into the future without making the optimization too slow.

I tried `N=30, dt=0.1` and `10, 0.05`. `N=30, dt=0.1` causes the vehicle to behave radically during sharp turns, 
and sometimes drive off the road, I think because the control it provides is too coarse. Whereas, `10, 0.05` causes 
the vehicle to slowly drift off the track, I think because when you can't see far enough into the future, you'll always
find yourself "surprised" and what's coming, but due to the fact the actuators are penalized by the cost function, there's no 
"incentive" for the optimizer to make sudden changes in the actuators.

# Polynomial Fitting and MPC Preprocessing

The waypoints are fit to 3rd order polynomials. And the waypoints are transformed to vehicle coordinate prior to fitting.
Thus for the states, `x`, `y`, and `psi` are always 0.

# Model Predictive Control with Latency

I deal with latency by simulating it by fixing the actuators to the previously returned values in the first two `dt` time period 
in the optimization cycle, so that in those two periods, the vehicle is assumed to always assume the last actuators values 
returned.
