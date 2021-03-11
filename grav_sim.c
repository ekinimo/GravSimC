#include <SDL2/SDL_events.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_timer.h>
#include <SDL2/SDL_video.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <SDL2/SDL.h>


#define SIM_DT (1e-12)
#define SIM_G  (10) 
#define WIDTH   (1200)
#define HEIGHT  (900)
#define ITER_UNTIL_UPDATE (10)


/*
TODO:
---- a Simulation state structure

---- collision detection is not physical.
------- Problem : static collision is unrealistic you should calculate the
-------           timepoint where two circles intersect
-------           and use that timepoint in dynamiccollision detection to calculate where
-------           and what speed the particle would have but that brings us to the next problem
------- Problem:  currently collisions are checked one by one and resolved in order. 
-------           that solves problems falsely. one should first record all collision and 
-------           then resolve them synchronously.

---- inelastic collisions
---- also two planets can crash and result in 3 new particles. add that to simulation

---- make it 3d

---- quadtree octotree would significantly improve the performance

---- calculate auxilary variables energy etc.

---- forget the planets that are gone rogue / what is the condition for rogue planet

---- proper draw circle / filled circle : midpoint algo

---- other integration schemes ?
--------Higher order symplectic
--------Runge-kutta
--------multi grid / gridless

---- charge ? spin ? 

*/

typedef struct Particle{
  size_t id;
  double mass;
  double radius;
  double x0,x1;
  double v0,v1;
  double a0,a1;
  uint8_t color[4];
  
}Particle_t;

typedef struct ParticleTree{
  Particle_t particle;
  Particle_t* left;
  Particle_t* right;
}ParticleTree_t;


static void print_particle(const Particle_t *p1) {
  printf("Particle ");
  printf("id :\t%zu\n",p1->id);
  printf("m/r:\t%f\t%f\t\n",p1->mass,p1->radius);
  printf("pos:\t%f\t%f\t\n",p1->x0,p1->x1);
  printf("vel:\t%f\t%f\t\n",p1->v0,p1->v1);
  printf("acc:\t%f\t%f\t\n",p1->a0,p1->a1);
  printf("\n-----------------------------------\n");
}

static void print_particle_array( Particle_t ** particles,size_t n) {
  printf("+++++++++++++++++++++PARTICLE_LIST++++++++++++++++++++++++\n");
  for(size_t i = 0 ; i < n ; i++)
    {
      print_particle(particles[i]);
    }
  printf("\n");
  
}

static uint8_t DoParticlesOverlap(const Particle_t* p1, const Particle_t* p2) {
  double d0 =  p1->x0 - p2->x0;
  double d1 =  p1->x1 - p2->x1;
  double r = p1->radius + p2->radius;
  return fabs(d0*d0 + d1*d1) <= r*r;
}

static double DistanceBetweenParticles(const Particle_t* p1, const Particle_t* p2) {
  double d0 =  p1->x0 - p2->x0;
  double d1 =  p1->x1 - p2->x1;
  return sqrt(d0*d0 + d1*d1);
}

static void StaticCollisionResolution(Particle_t* p1, Particle_t* p2) {
  double distance = DistanceBetweenParticles(p1,p2);
  double overlap = 0.5*(distance - p1->radius - p2->radius);
  // static collision resolution
  p1->x0 -= overlap *(p1->x0 - p2->x0) / distance;
  p1->x1 -= overlap *(p1->x1 - p2->x1) / distance;
  
  p2->x0 += overlap *(p1->x0 - p2->x0) / distance;
  p2->x1 += overlap *(p1->x1 - p2->x1) / distance;

}
static void DynamicCollisionResolution(Particle_t* p1, Particle_t* p2) {

  double distance = DistanceBetweenParticles(p1,p2);

  //normal vector
  double n0 =  (p1->x0 - p2->x0)/distance;
  double n1 =  (p1->x1 - p2->x1)/distance;

  //tangent vector
  double t0 = -n1;
  double t1 = n0;

  double dtan0 = p1->v0 * t0 + p1->v1 * t1;
  double dtan1 = p2->v0 * t0 + p2->v1 * t1;
  double dnorm0 = p1->v0 * n0 + p1->v1 * n1;
  double dnorm1 = p2->v0 * n0 + p2->v1 * n1;

  double m1 = (dnorm0 * (p1->mass - p2->mass) + 2 * p2->mass *dnorm1)/(p1->mass+p2->mass);
  double m2 = (dnorm1 * (p1->mass - p2->mass) + 2 * p1->mass *dnorm0)/(p1->mass+p2->mass);

  p1->v0 = t0 * dtan0 + m1* n0  ;
  p1->v1 = t1 * dtan0 + m1* n1 ;
  p2->v0 = t0 * dtan1 + m2* n0;
  p2->v1 = t1 * dtan1 + m2* n1;
  
}

static void collide(Particle_t *p1,  Particle_t *p2) {

  if(p1->id != p2->id)
    {
      if (DoParticlesOverlap(p1,p2))
	{
	  printf("DETECTED COLLISSIION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	  StaticCollisionResolution(p1,p2);
	  DynamicCollisionResolution(p1,p2);
	  
	}
    }
  
}

static void collision_naive(Particle_t *p1,
			     Particle_t  ** particles,
			    size_t n) {

  for(size_t i = 0 ; i < p1->id ; i++)
    {
      Particle_t* p2 = particles[i]; 
      collide(p1,p2);
    }
}

static void calculate_accelaration_naive(Particle_t *p1,
                                            Particle_t *  * particles,
					   size_t n) {
  
  p1->a0 = 0;
  p1->a1 = 0;
  
  for(size_t i = 0 ; i < n ; i++)
    {
      
       Particle_t* p2 = particles[i]; 
       if(p2->id != p1->id){
	 double dist = DistanceBetweenParticles(p1,p2);
	 double temp= SIM_G * p2->mass / (dist * dist * dist);
	 p1->a0 += -(p1->x0 - p2->x0)*temp;
	 p1->a1 += -(p1->x1 - p2->x1)*temp;
       }
    }
}

static void leapfrog_naive(Particle_t *p1,
		      Particle_t *  * particles,
		     size_t n) {

  calculate_accelaration_naive(p1,particles,n);
  p1->v0 += p1->a0 + SIM_DT /2 ;
  p1->v1 += p1->a1 + SIM_DT /2 ;

  p1->x0 += p1->v0 + SIM_DT  ;
  p1->x1 += p1->v1 + SIM_DT  ;

  calculate_accelaration_naive(p1,particles,n);
  p1->v0 += p1->a0 + SIM_DT /2 ;
  p1->v1 += p1->a1 + SIM_DT /2 ;
  
}

static void evolve(  Particle_t **particles, size_t n) {
  for(size_t i = 0 ; i < n ; i++)
    {
      collision_naive(particles[i],particles,n);
      leapfrog_naive(particles[i],particles,n);
      
    }  
}

Particle_t *new_particle(double radius,double mass,
			 double x0,double x1,
			 double v0,double v1) {
  Particle_t *ret = calloc(1,sizeof(Particle_t));
  if(ret ==NULL){printf("couldnt allocate particle exiting \n"); exit(0);}
  static size_t id ;
  ret->id = id++ ;
  ret->mass = mass; ret->radius = radius;
      ret->x0 = x0;  ret->x1 = x1 ;
      ret->v0 = v0;    ret->v1= v1 ;
      ret->a0=0;      ret->a1=0;
      ret->color[0] = rand();
      ret->color[1] = rand();
      ret->color[2] = rand();
      ret->color[3] = 255;
      
  return ret;
}




static void draw_reg_polygon(SDL_Renderer* renderer , float x, float y,float r,int n) {

  const float incr =  (2 * M_PI) /((float) n) ;

  float p1 = 0 ;
  float p2 = 0;

  float p1_old = x ;
  float p2_old = y + r;

  for(float i = 1 ; i<= n ; i++ ){
    p1 = x + r*sin( i* incr);
    p2 = y + r*cos( i* incr);
    SDL_RenderDrawLineF(renderer,p1,p2,p1_old,p2_old);
    p1_old = p1;
    p2_old = p2;
  }
}

static void draw_particles(SDL_Renderer *renderer , Particle_t **particles, size_t n) {

  for(size_t i = 0 ; i < n ; i++)
    {
      
      float x = particles[i]->x0;
      float y = particles[i]->x1;
      float r = particles[i]->radius;
      uint8_t* col = particles[i]->color;
      SDL_SetRenderDrawColor(renderer,col[0],col[1],col[2],col[3]);
      draw_reg_polygon(renderer,
		       x,
		       y,
		       r,20);
      //printf("\t\n\n%f\t%f\t%f\n",x,y,r);
    }
  
}

int main() {

  size_t num_of_particles =  3;
  Particle_t ** particles = calloc(num_of_particles,sizeof(Particle_t*));
  particles[0] = new_particle(10,50,WIDTH/2,HEIGHT/2,0,0);
  particles[1] = new_particle(5,1,WIDTH/2-100,HEIGHT/2,0,0);
  particles[2] = new_particle(5,1,WIDTH/2+100,HEIGHT/2,0,0);
  /*for(size_t i = 1 ;i<num_of_particles;i++){
    double mass = abs(rand()) % 100;
    double rad = abs(rand()) % 30;

    double x0 = rand() % HEIGHT;
    double x1 = rand() % WIDTH;
    double v0 = (((float)(rand() % 10)) - 5)/2;
    double v1 = (((float)(rand() % 10)) - 5)/2;

    particles[i] = new_particle(rad,mass,x0,x1,v0,v1);
  }
  */
  if (SDL_Init(SDL_INIT_EVERYTHING)!=0)
    {
      printf("couldn init sdl exiting");exit(0);
    }

  SDL_Window* window;
  SDL_Renderer* renderer;
  if (SDL_CreateWindowAndRenderer(WIDTH,HEIGHT,SDL_WINDOW_MOUSE_CAPTURE,&window,&renderer)!=0)
    {
      printf("couldn init sdl exiting");exit(0);
    }

  int running = 1;
  SDL_Event event;
  SDL_SetRenderDrawColor(renderer,0,0,0,255); 
  SDL_RenderClear(renderer);
  while(running)
    {
      while(SDL_PollEvent(&event))
	{
	  switch(event.type)
	    {
	    case SDL_KEYDOWN:
	      if (event.key.keysym.sym == SDLK_q){running=0;}
	      break;
	    case SDL_QUIT:
	      running = 0;
	      break;
	    default:
	      break;
	    }
	}
      //for(int i = 0 ; i < ITER_UNTIL_UPDATE; i++)
      //{
	  evolve(particles,num_of_particles);
      //}
      SDL_SetRenderDrawColor(renderer,0,0,0,255);
      SDL_RenderClear(renderer);
      SDL_SetRenderDrawColor(renderer,255,255,255,255);

      SDL_Delay(50);
      draw_particles(renderer,particles,num_of_particles);
      //print_particle_array(particles,2);
      SDL_RenderPresent(renderer);
      
    }
  
}
