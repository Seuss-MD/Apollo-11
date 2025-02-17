/*************************************************************
* 1. Name:
*      -Sam Evans, Savannah Harvey-
* 2. Assignment Name:
*      Practice 02: Physics simulator
* 3. Assignment Description:
*      Compute how the Apollo lander will move across the screen
* 4. What was the hardest part? Be as specific as possible.
*      Understanding the physics for the assignment
* 5. How long did it take for you to complete the assignment?
*      About 4 hours
**************************************************************/

#include <iostream>  // for CIN and COUT
#include <cmath>     // for pi, sin(), and cos()
using namespace std;

#define WEIGHT   15103.000   // Weight in KG
#define GRAVITY     -1.625   // Vertical acceleration due to gravity, in m/s^2
#define THRUST   45000.000   // Thrust of main engine, in Newtons (kg m/s^2)

/***************************************************
* COMPUTE DISTANCE
* Apply inertia to compute a new position using the distance equation.
* The equation is:
*     s = s + v t + 1/2 a t^2
* INPUT
*     s : original position, in meters
*     v : velocity, in meters/second
*     a : acceleration, in meters/second^2
*     t : time, in seconds
* OUTPUT
*     s : new position, in meters
**************************************************/
double computeDistance(double s, double v, double a, double t)
{
   // get updated position
   s = s + (v * t) + (1.0 / 2.0) * a* t* t;

   return s;
}
 

/**************************************************
* COMPUTE ACCELERATION
* Find the acceleration given a thrust and mass.
* This will be done using Newton's second law of motion:
*     f = m * a
* INPUT
*     f : force, in Newtons (kg * m / s^2)
*     m : mass, in kilograms
* OUTPUT
*     a : acceleration, in meters/second^2
***************************************************/
double computeAcceleration(double f, double m)
{
   double acceleration;

   // computes acceleration using a = f / m
   acceleration = f / m;

   return acceleration;
}


/***********************************************
* COMPUTE VELOCITY
* Starting with a given velocity, find the new
* velocity once acceleration is applied. This is
* called the Kinematics equation. The
* equation is:
*     v = v + a t
* INPUT
*     v : velocity, in meters/second
*     a : acceleration, in meters/second^2
*     t : time, in seconds
* OUTPUT
*     v : new velocity, in meters/second
***********************************************/
double computeVelocity (double velocity, double acceleration, double time)
{

   //computes the velocity
   velocity = velocity + acceleration * time;

   return velocity;
}


/***********************************************
* COMPUTE VERTICAL COMPONENT
* Find the vertical component of a velocity or acceleration.
* The equation is:
*     cos(a) = y / total
* This can be expressed graphically:
*      x
*    +-----
*    |   /
*  y |  / total
*    |a/
*    |/
* INPUT
*     a : angle, in radians
*     total : total velocity or acceleration
* OUTPUT
*     y : the vertical component of the total
***********************************************/
double computeVertical(double angle, double total)
{
   double verticalY;

   // compute y, y = cos(a) * total
   verticalY = cos(angle) * total;

   return verticalY;
}

/***********************************************
* COMPUTE HORIZONTAL COMPONENT
* Find the horizontal component of a velocity or acceleration.
* The equation is:
*     sin(a) = x / total
* This can be expressed graphically:
*      x
*    +-----
*    |   /
*  y |  / total
*    |a/
*    |/
* INPUT
*     a : angle, in radians
*     total : total velocity or acceleration
* OUTPUT
*     x : the vertical component of the total
***********************************************/
double computeHorizontal(double a, double total)
{
   double x;

   // x = sin(a) * total
   x = sin(a) * total;

   return x;
}

/************************************************
* COMPUTE TOTAL COMPONENT
* Given the horizontal and vertical components of
* something (velocity or acceleration), determine
* the total component. To do this, use the Pythagorean Theorem:
*    x^2 + y^2 = t^2
* where:
*      x
*    +-----
*    |   /
*  y |  / total
*    | /
*    |/
* INPUT
*    x : horizontal component
*    y : vertical component
* OUTPUT
*    total : total component
***********************************************/
double computeTotal(double x, double y)
{
   double total;

   // total = square root(x^2 + y^2)
   total = sqrt((x * x) + (y * y));

   return total;
}


/*************************************************
* RADIANS FROM DEGEES
* Convert degrees to radians:
*     radians / 2pi = degrees / 360
* INPUT
*     d : degrees from 0 to 360
* OUTPUT
*     r : radians from 0 to 2pi
**************************************************/
double getRadian(double d)
{
   // calculate pi
   double pi = 2.0 * asin(1.0);

   //convert degrees to radians
   double r;
   r = d / 360.0 * 2.0 * pi;

   return r;
}

/**************************************************
* PROMPT
* A generic function to prompt the user for a double
* INPUT
*      message : the message to display to the user
* OUTPUT
*      response : the user's response
***************************************************/
double prompt(string message)
{
   //display message
   cout << message;

   //get response
   double response;
   cin  >> response;

   return response;
}

/****************************************************************
 * MAIN
 * Prompt for input, compute new position, and display output
 ****************************************************************/
int main()
{
   // Prompt for input and variables to be computed
   double dx = prompt("What is your horizontal velocity (m/s)? ");
   double dy = prompt("What is your vertical velocity (m/s)? ");
   double y = prompt("What is your altitude (m)? ");
   double x = prompt("What is your position (m)? ");
   double aDegrees = prompt("What is the angle of the LM where 0 is up (degrees)? ");
   double t = prompt("What is the time interval (s)? ");
   double aRadians;            // Angle in radians
   double accelerationThrust;  // Acceleration due to thrust 
   double ddxThrust;           // Horizontal acceleration due to thrust
   double ddyThrust;           // Vertical acceleration due to thrust
   double ddx;                 // Total horizontal acceleration
   double ddy;                 // Total vertical acceleration
   double v;                   // Total velocity

   // get accelerationThrust
   accelerationThrust = computeAcceleration(THRUST, WEIGHT);

   // Go through the simulator five times
   for (int i = 0; i < 5; i++)
   {
      // get radians
      aRadians = getRadian(aDegrees);
      
      // get the total acceleration of the horizontal and vertical
      ddx = computeHorizontal(aRadians, accelerationThrust);
      ddy = computeVertical(aRadians, accelerationThrust) + GRAVITY;
      

      //get the new position
      y = computeDistance(y, dy, ddy, t);
      x = computeDistance(x, dx, ddx, t);
      
      // get the new velocity
      dx = computeVelocity(dx, ddx, t);
      dy = computeVelocity(dy, ddy, t);
      v  = computeTotal(dx, dy);

      // Output
      cout.setf(ios::fixed | ios::showpoint);
      cout.precision(2);
      cout << "\tNew position:   (" << x << ", " << y << ")m\n";
      cout << "\tNew velocity:   (" << dx << ", " << dy << ")m/s\n";
      cout << "\tTotal velocity:  " << v << "m/s\n\n";
   }

   return 0;
}