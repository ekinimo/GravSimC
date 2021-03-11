# Gravity Simulator #
Hey this is an ongoing project. I aim to implement Barnes Hut algorithm in the end with elastic/inelastic collisions. We will se how it turns out

# TODO #
  * a Simulation state structure
  *  collision detection is not physical.
    *   Problem : static collision is unrealistic you should calculate the timepoint where two circles intersect and use that timepoint in dynamiccollision detection to calculate where and what speed the particle would have but that brings us to the next problem
	
    *   Problem:  currently collisions are checked one by one and resolved in order. That solves problems falsely. one should first record all collision and then resolve them synchronously.


  * inelastic collisions

  *  also two planets can crash and result in 3 new particles. add that to simulation


  *  make it 3d


  *  quadtree octotree would significantly improve the performance


  *  calculate auxilary variables energy etc.


  *  forget the planets that are gone rogue / what is the condition for rogue planet


  *  proper draw circle / filled circle : midpoint algo


  *  other integration schemes ? 
      *  Higher order symplectic
      *  Runge-kutta
      *  multi grid / gridless


  * charge ? spin ? 
