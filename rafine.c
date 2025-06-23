#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define FILENUMBER 1000     // Number appended to end of file
#define DATAFILE "data"              // Filenames
#define HINGEFILE "hingedata"
#define RODFILE "roddata"
#define MIDPOINTFILE "midpointdata"
#define STICKFILE "stickdata"

#define SHEARDISPL 0.01L     // Only one of these should be nonzero
#define STRETCHDISPL 0.0L

#define WX 1.0L // Size of domain in x direction (m)
#define WY 1.0L // Size of domain in y direction (m)
#define L 0.2L  // Length of rods (m) (Must have L<min(WX/2,WY/2) => only
               // one intersection between rods)
#define N 500 // Number of rods (1)
#define EA 1.0L // Normal stiffness of rods (N)
#define EI 0.00000004L // Bending stiffness (Nm^2)
#define GI 0.0000000001L // Shear stiffness (Nm)
#define SEED 105 // Seed for random number generator
#define CASENO 10 // Number of runs with the same parameters
#define SMALL 1.0e-12L

// **************************
// Data structures
// **************************

struct rods  /* A chain containing all rods */
{
	long double x;    // Central point, x coordinate
	long double y;    // Central point, x coordinate
	long double teta; // Angle with horizontal, counter clockwise
	long double tt;   // tan(teta)
	long double sit;  // sin(teta)
	long double cot;  // cos(teta)
	long double x0;   // End points
	long double y0;
	long double x1;
	long double y1;
	struct rods *next;
};

struct hinges /* A chain containing all hinges */
{
	long int ser;   // Serial number of hinge (0 is hinge is on edge!)
	long double x0;  // Original x coordinate
	long double y0;  // Original y coordinate
	long double x;   // Actual x coordinate
	long double y;   // Actual y coordinate
	struct rods *rod1;  // One of the rods intersecting at this hinge
	struct rods *rod2;  // Other rod intersecting at this hinge
	struct sticks *st1;  // Pointer to neighbouring sticks (2,3 or 4),
	struct sticks *st2;  // st1 and st2 are on rod1, st3 and st4 are on
	struct sticks *st3;  // rod 2. If one of these is 0, no neighbouring
	struct sticks *st4;  // stick there.
	struct hinges *n1;  // Pointer to neighbouring hinges (2,3 or 4),
	struct hinges *n2;  // n1 and n2 are on rod1, n3 and n4 are on rod2.
	struct hinges *n3;  // If one of these is 0, no neighbouring hinge
	struct hinges *n4;  // there.
	struct midpoints *mp1;  // Pointer to neighbouring midpoints (2,3 or 4),
	struct midpoints *mp2;  // mp1 and mp2 are on rod1, mp3 and mp4 are on
	struct midpoints *mp3;  // rod2. If one of these is 0, no neighbouring
	struct midpoints *mp4;  // midpoint there.
	struct hinges *next;
};

struct midpoints /* A chain containing all midpoints of sticks */
{
	long int ser;       // Serial number of midpoint
	long double x0;     // Original x coordinate
	long double y0;     // Original y coordinate
	long double x;      // Actual x coordinate
	long double y;      // Actual y coordinate
	int n1x;            // These indicate if hinge n1 or n2 is out of domain:
	int n1y;            //     0 for inside, +1 for hinge out in positive x or y
	int n2x;            //     direction, -1 for hinge out in negative x or y
	int n2y;            //     direction.
	struct sticks *st;  // Pointer to stick containing this midpoint
	struct hinges *n1;  // Pointer to neighbouring hinges
	struct hinges *n2;
	struct midpoints *next;
};

struct sticks /* A chain containing all rod segments between hinges */
{
	struct hinges *n1;      // One end
	struct hinges *n2;      // Other end
	struct midpoints *mp;  // Midpoint
	struct rods *rod;       // Rod containing this stick
	long double len;             // Length of stick in original state
	struct sticks *next;
};

// ****************************
// Global variables
// ****************************

struct rods *rodstart;
struct hinges *hingestart;
struct midpoints *midpointstart;
struct sticks *stickstart;
FILE *data, *roddata, *hingedata, *midpointdata, *stickdata;
long double *bv;
long double cconst;
long double PI=acosl(-1.0L);


long int kspr,ksprMAX,ksprSTP;
long double *Aspr; 
long int *Aspri,*Asprj;  // SPARSE i,j indexek

// ********************
// Subprograms
// ********************

int mkAsp(long int SS, long int k, long double s)
//fills the Aspr matrix
{
	long int err;
	long int l; long int ks;
	long int i; long int j;
	j= (k%SS)+1;
	i= (k/SS)+1;
	err=0;
// ha meg nincs Aspr, akkor le kell foglalni
	if(ksprMAX==0)
	{
			ksprMAX=ksprMAX+ksprSTP;
			Aspr=(long double *)calloc(ksprMAX,sizeof(long double));
			Aspri=(long int *)calloc(ksprMAX,sizeof(long int));
			Asprj=(long int *)calloc(ksprMAX,sizeof(long int));
	}
//vegig kell nezni 1-kspr-ig, hogy van-e mar i,j index
//ha van i.j index, akkor az lesz a ks
	ks=0;
	for (l=1;l<=kspr;l++)
	{
		if(Aspri[l]==i && Asprj[l]==j)
		{ks=l;
		}
	}
//ha nincs meg i,j index, akkor a kspr-t növeljük és az lesz a ks
	if(ks==0)
	{
		kspr=kspr+1; ks=kspr;
//ha még nincs annyi lefoglalva, akkor növelni a lefoglalast
		if(ksprMAX==kspr)
		{
			ksprMAX=ksprMAX+ksprSTP;
			Aspr=(long double *)realloc(Aspr,ksprMAX*sizeof(long double));
			Aspri=(long int *)realloc(Aspri,ksprMAX*sizeof(long int));
			Asprj=(long int *)realloc(Asprj,ksprMAX*sizeof(long int));
			if(Aspr == NULL){  printf("OOM\n");}
			for(l=ksprMAX-ksprSTP;l<ksprMAX;l++)
			{
			Aspr[l]=0.0L; //es inicializaljuk is
			Aspri[l]=0;
			Asprj[l]=0;
			}
		}
		Aspri[ks]=i;
		Asprj[ks]=j;
	}
	Aspr[ks]=Aspr[ks]+s;
	return err;
}

int setuprods()
// Generates chain of rods
{
	int i;
	struct rods *p;
	rodstart=(struct rods *)calloc(1,sizeof(struct rods));
	p=rodstart;
	for (i=0;i<N-1;i++)
	{
		(*p).x=drand48()*WX; // random, 0-WX
		(*p).y=drand48()*WY; // random, 0-WY
		(*p).teta=drand48()*PI-PI/2.0; // random, (-pi/2)-(+pi/2)
		(*p).tt=tanl((*p).teta);
		(*p).sit=sinl((*p).teta);
		(*p).cot=cosl((*p).teta);
		(*p).x0=(*p).x-L*cosl((*p).teta)/2.0L;
		(*p).y0=(*p).y-L*sinl((*p).teta)/2.0L;
		(*p).x1=(*p).x+L*cosl((*p).teta)/2.0L;
		(*p).y1=(*p).y+L*sinl((*p).teta)/2.0L;
		(*p).next=(struct rods *)calloc(1,sizeof(struct rods));
		p=(*p).next;
	}
	(*p).x=drand48()*WX; // random, 0-WX
	(*p).y=drand48()*WY; // random, 0-WY
	(*p).teta=drand48()*PI-PI/2.0L; // random, (-pi/2)-(+pi/2)
	(*p).tt=tanl((*p).teta);
	(*p).sit=sinl((*p).teta);
	(*p).cot=cosl((*p).teta);
	(*p).x0=(*p).x-L*cosl((*p).teta)/2.0L;
	(*p).y0=(*p).y-L*sinl((*p).teta)/2.0L;
	(*p).x1=(*p).x+L*cosl((*p).teta)/2.0L;
	(*p).y1=(*p).y+L*sinl((*p).teta)/2.0L;
	(*p).next=0;
	return 0;
}

// **************************

int intersect(long double x0, long double y0, long double x1, long double y1, long double t0, long double t1, long double x0s, long double y0s, long double x0e, long double y0e, long double x1s, long double y1s, long double x1e, long double y1e, long double *x, long double *y)
// Finds intersection of two line segments
{
	int err;
	long double dum;
	// Assume there is no intersection 
	err=0;
	if (fabsl(t0-t1)>SMALL)
	{
		// Intersection of lines
		*x=(x0*t0-x1*t1+y1-y0)/(t0-t1);
		*y=((x0-x1)*t0*t1+y1*t0-y0*t1)/(t0-t1);
		// See if intersection is on both segments
		if (x0s>x0e)
		{
			dum=x0s;
			x0s=x0e;
			x0e=dum;
		}
		if (x1s>x1e)
		{
			dum=x1s;
			x1s=x1e;
			x1e=dum;
		}
		if (y0s>y0e)
		{
			dum=y0s;
			y0s=y0e;
			y0e=dum;
		}
		if (y1s>y1e)
		{
			dum=y1s;
			y1s=y1e;
			y1e=dum;
		}
		if ((*x<=x0e+SMALL)&&(*x+SMALL>=x0s))
		{
			if ((*x<=x1e+SMALL)&&(*x+SMALL>=x1s))
			{
				if ((*y<=y0e+SMALL)&&(*y+SMALL>=y0s))
				{
					if ((*y<=y1e+SMALL)&&(*y+SMALL>=y1s))
					{
						// Segments intersect
						err=1;
					}
				}
			}
		}
	}
	if (err==1)
	{
		while (*x<0.0) *x=*x+WX;
		while (*x>WX) *x=*x-WX;
		while (*y<0.0) *y=*y+WY;
		while (*y>WY) *y=*y-WY;
	}
	return err;
}

// **************************

int findhinges()
// Finds all intersections between rods and perimeter. Output is number of hinges found.
{
	struct rods *p,*q;
	struct hinges *g,*h;
	int err,hno;
	long double outx,outy,outxo,outyo;
	long double *x,*y;
	x=(long double *)calloc(1,sizeof(long double));
	y=(long double *)calloc(1,sizeof(long double));
	p=rodstart;
	hingestart=(struct hinges *)calloc(1,sizeof(struct hinges));
	h=hingestart;
	g=h;
	hno=0;
	while (p)
	{
		// Does this rod hang over the edge?
		// If yes, hinges are introduced on opposite edges
		outx=0.0L;
		outy=0.0L;
		if (((*p).x0<0.0L)||((*p).x1<0.0L))
		{
			// Hangs out on -x side
			outx=-WX;
		}
		if (((*p).x0>WX)||((*p).x1>WX))
		{
			// Hangs out on +x side
			outx=WX;
		}
		if (((*p).y0<0.0L)||((*p).y1<0.0L))
		{
			// Hangs out on -y side
			outy=-WY;
		}
		if (((*p).y0>WY)||((*p).y1>WY))
		{
			// Hangs out on +y side
			outy=WY;
		}
		// Cycle along all subsequent rods to look for intersections
		q=(*p).next;
		while (q)
		{
			// See if this other rod hangs out
			outxo=0.0L;
			outyo=0.0L;
			if (((*q).x0<0.0L)||((*q).x1<0.0L))
			{
				// Hangs out on -x side
				outxo=-WX;
			}
			if (((*q).x0>WX)||((*q).x1>WX))
			{
				// Hangs out on +x side
				outxo=WX;
			}
			if (((*q).y0<0.0L)||((*q).y1<0.0L))
			{
				// Hangs out on -y side
				outyo=-WY;
			}
			if (((*q).y0>WY)||((*q).y1>WY))
			{
				// Hangs out on +y side
				outyo=WY;
			}
			// Go from case to case: find intersection
			if ((fabsl(outx)<SMALL)&&(fabsl(outy)<SMALL))
			{
				// first rod is inside totally
				err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
				// err: 0 if no intersection, 1 if there is,
				// *x,*y has the intersection.
				if (err==1)
				{
					// New hinge
					hno++;
					(*h).ser=hno;
					(*h).x0=*x;
					(*h).y0=*y;
					(*h).x=0.0L;
					(*h).y=0.0L;
					(*h).rod1=p;
					(*h).rod2=q;
					(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
					h=(*h).next;
					if (!(h==(*hingestart).next)) g=(*g).next;
				}
				else
				{
					// Is there a hinge otherwise? Push second rod back on both sides
					err=intersect((*p).x,(*p).y,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
					if (err==1)
					{
						// New hinge
						hno++;
						(*h).ser=hno;
						(*h).x0=*x;
						(*h).y0=*y;
						(*h).x=0.0L;
						(*h).y=0.0L;
						(*h).rod1=p;
						(*h).rod2=q;
						(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
						h=(*h).next;
						if (!(h==(*hingestart).next)) g=(*g).next;
					}
					else
					{
						// Only remaining possibility for hinge is that second rod
						// hangs out on two sides, and needs to be pushed back on
						// one side only.
						if ((fabsl(outxo)>SMALL)&&(fabsl(outyo)>SMALL))
						{
							// Second rod hangs out on two sides
							// Pull it back along x:
							err=intersect((*p).x,(*p).y,(*q).x-outxo,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0-outxo,(*q).y0,(*q).x1-outxo,(*q).y1,x,y);
							if (err==1)
							{
								// New hinge
								hno++;
								(*h).ser=hno;
								(*h).x0=*x;
								(*h).y0=*y;
								(*h).x=0.0L;
								(*h).y=0.0L;
								(*h).rod1=p;
								(*h).rod2=q;
								(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
								h=(*h).next;
								if (!(h==(*hingestart).next)) g=(*g).next;
							}
							else
							{
							// Pull it back along y:
								err=intersect((*p).x,(*p).y,(*q).x,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0-outyo,(*q).x1,(*q).y1-outyo,x,y);
								if (err==1)
								{
									// New hinge
									hno++;
									(*h).ser=hno;
									(*h).x0=*x;
									(*h).y0=*y;
									(*h).x=0.0L;
									(*h).y=0.0L;
									(*h).rod1=p;
									(*h).rod2=q;
									(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
									h=(*h).next;
									if (!(h==(*hingestart).next)) g=(*g).next;
								}
							}
						}

					}
				}
			}
			else
			{
				if ((fabsl(outx)>SMALL)&&(fabsl(outy)>SMALL))
				{
					// first rod hangs out on two sides.
					// See what the second rod does.
					// Is it totally inside?
					if ((fabsl(outxo)<SMALL)&&(fabsl(outyo)<SMALL))
					{
						// Second rod is totally inside
						err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
						if (err==1)
						{
							// New hinge
							hno++;
							(*h).ser=hno;
							(*h).x0=*x;
							(*h).y0=*y;
							(*h).x=0.0L;
							(*h).y=0.0L;
							(*h).rod1=p;
							(*h).rod2=q;
							(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
							h=(*h).next;
							if (!(h==(*hingestart).next)) g=(*g).next;
						}
						else
						{
							err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
							if (err==1)
							{
								// New hinge
								hno++;
								(*h).ser=hno;
								(*h).x0=*x;
								(*h).y0=*y;
								(*h).x=0.0L;
								(*h).y=0.0L;
								(*h).rod1=p;
								(*h).rod2=q;
								(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
								h=(*h).next;
								if (!(h==(*hingestart).next)) g=(*g).next;
							}
							else
							{
								err=intersect((*p).x-outx,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0,(*p).x1-outx,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
								if (err==1)
								{
									// New hinge
									hno++;
									(*h).ser=hno;
									(*h).x0=*x;
									(*h).y0=*y;
									(*h).x=0.0L;
									(*h).y=0.0L;
									(*h).rod1=p;
									(*h).rod2=q;
									(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
									h=(*h).next;
									if (!(h==(*hingestart).next)) g=(*g).next;
								}
								else
								{
									err=intersect((*p).x,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0-outy,(*p).x1,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
									if (err==1)
									{
										// New hinge
										hno++;
										(*h).ser=hno;
										(*h).x0=*x;
										(*h).y0=*y;
										(*h).x=0.0L;
										(*h).y=0.0L;
										(*h).rod1=p;
										(*h).rod2=q;
										(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
										h=(*h).next;
										if (!(h==(*hingestart).next)) g=(*g).next;
									}
								}
							}
						}
					}
					else
					{
						if ((fabsl(outxo*outyo)<SMALL)&&((fabsl(outx-outxo)<SMALL)||(fabsl(outy-outyo)<SMALL)))
						{
							// Second rod is out on one side, in the same direction as the first
							err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
							if (err==1)
							{
								// New hinge
								hno++;
								(*h).ser=hno;
								(*h).x0=*x;
								(*h).y0=*y;
								(*h).x=0.0L;
								(*h).y=0.0L;
								(*h).rod1=p;
								(*h).rod2=q;
								(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
								h=(*h).next;
								if (!(h==(*hingestart).next)) g=(*g).next;
							}
							else
							{
								err=intersect((*p).x-outx+outxo,(*p).y-outy+outyo,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx+outxo,(*p).y0-outy+outyo,(*p).x1-outx+outxo,(*p).y1-outy+outyo,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
								if (err==1)
								{
									// New hinge
									hno++;
									(*h).ser=hno;
									(*h).x0=*x;
									(*h).y0=*y;
									(*h).x=0.0L;
									(*h).y=0.0L;
									(*h).rod1=p;
									(*h).rod2=q;
									(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
									h=(*h).next;
									if (!(h==(*hingestart).next)) g=(*g).next;
								}
							}
						}
						else
						{
							if ((fabsl(outxo*outyo)<SMALL)&&((fabsl(outx+outxo)<SMALL)||(fabsl(outy+outyo)<SMALL)))
							{
								// Second rod is out on a side where first rod is not out
								err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
								if (err==1)
								{
									// New hinge
									hno++;
									(*h).ser=hno;
									(*h).x0=*x;
									(*h).y0=*y;
									(*h).x=0.0L;
									(*h).y=0.0L;
									(*h).rod1=p;
									(*h).rod2=q;
									(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
									h=(*h).next;
									if (!(h==(*hingestart).next)) g=(*g).next;
								}
								else
								{
									err=intersect((*p).x,(*p).y,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
									if (err==1)
									{
										// New hinge
										hno++;
										(*h).ser=hno;
										(*h).x0=*x;
										(*h).y0=*y;
										(*h).x=0.0L;
										(*h).y=0.0L;
										(*h).rod1=p;
										(*h).rod2=q;
										(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
										h=(*h).next;
										if (!(h==(*hingestart).next)) g=(*g).next;
									}
								}
							}
							else
							{
								if ((fabsl(outxo*outyo)>SMALL)&&((fabsl(outx-outxo)<SMALL)&&(fabsl(outy-outyo)<SMALL)))
								{
									// Second rod is out on the same two sides as the first rod
									err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
									if (err==1)
									{
										// New hinge
										hno++;
										(*h).ser=hno;
										(*h).x0=*x;
										(*h).y0=*y;
										(*h).x=0.0L;
										(*h).y=0.0L;
										(*h).rod1=p;
										(*h).rod2=q;
										(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
										h=(*h).next;
										if (!(h==(*hingestart).next)) g=(*g).next;
									}
								}
								else
								{
									if ((fabsl(outxo*outyo)>SMALL)&&(((fabsl(outx-outxo)<SMALL)&&(fabsl(outy-outyo)>SMALL))||((fabsl(outx-outxo)>SMALL)&&(fabsl(outy-outyo)<SMALL))))
									{
										// Second rod is out on two sides, one is the same where first rod is out, the other is opposite
										err=intersect((*p).x-0.5L*(outx-outxo),(*p).y-0.5L*(outy-outyo),(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-0.5L*(outx-outxo),(*p).y0-0.5L*(outy-outyo),(*p).x1-0.5L*(outx-outxo),(*p).y1-0.5L*(outy-outyo),(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
										if (err==1)
										{
											// New hinge
											hno++;
											(*h).ser=hno;
											(*h).x0=*x;
											(*h).y0=*y;
											(*h).x=0.0L;
											(*h).y=0.0L;
											(*h).rod1=p;
											(*h).rod2=q;
											(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
											h=(*h).next;
											if (!(h==(*hingestart).next)) g=(*g).next;
										}
									}
									else
									{
										if ( (fabsl(outxo*outyo)>SMALL)&&((fabsl(outx+outxo)<SMALL)&&(fabsl(outy+outyo)<SMALL)))
										{
											// Second rod is out on two sides, both are opposite to sides where first rod is out
											err=intersect((*p).x-outx,(*p).y,(*q).x,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0,(*p).x1-outx,(*p).y1,(*q).x0,(*q).y0-outyo,(*q).x1,(*q).y1-outyo,x,y);
											if (err==1)
											{
												// New hinge
												hno++;
												(*h).ser=hno;
												(*h).x0=*x;
												(*h).y0=*y;
												(*h).x=0.0L;
												(*h).y=0.0L;
												(*h).rod1=p;
												(*h).rod2=q;
												(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
												h=(*h).next;
												if (!(h==(*hingestart).next)) g=(*g).next;
											}
											else
											{
												err=intersect((*p).x-outx,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0,(*p).x1-outx,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
												if (err==1)
												{
													// New hinge
													hno++;
													(*h).ser=hno;
													(*h).x0=*x;
													(*h).y0=*y;
													(*h).x=0.0L;
													(*h).y=0.0L;
													(*h).rod1=p;
													(*h).rod2=q;
													(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
													h=(*h).next;
													if (!(h==(*hingestart).next)) g=(*g).next;
												}
												else
												{
													err=intersect((*p).x,(*p).y,(*q).x,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0-outyo,(*q).x1,(*q).y1-outyo,x,y);
													if (err==1)
													{
														// New hinge
														hno++;
														(*h).ser=hno;
														(*h).x0=*x;
														(*h).y0=*y;
														(*h).x=0.0L;
														(*h).y=0.0L;
														(*h).rod1=p;
														(*h).rod2=q;
														(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
														h=(*h).next;
														if (!(h==(*hingestart).next)) g=(*g).next;
													}
												}
											}
										}	
									}
								}
							}
						}
					}
				}
				else
				{
					// first rod hangs out on one side
					// See what the second rod does.
					// Is it completely inside?
					if ((fabsl(outxo)<SMALL)&&(fabsl(outyo)<SMALL))
					{
						// Second rod is inside
						err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
						if (err==1)
						{
							// New hinge
							hno++;
							(*h).ser=hno;
							(*h).x0=*x;
							(*h).y0=*y;
							(*h).x=0.0L;
							(*h).y=0.0L;
							(*h).rod1=p;
							(*h).rod2=q;
							(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
							h=(*h).next;
							if (!(h==(*hingestart).next)) g=(*g).next;
						}
						else
						{
							err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
							if (err==1)
							{
								// New hinge
								hno++;
								(*h).ser=hno;
								(*h).x0=*x;
								(*h).y0=*y;
								(*h).x=0.0L;
								(*h).y=0.0L;
								(*h).rod1=p;
								(*h).rod2=q;
								(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
								h=(*h).next;
								if (!(h==(*hingestart).next)) g=(*g).next;
							}
						}
					}
					else
					{
						if ((fabsl(outx-outxo)<SMALL)&&(fabsl(outy-outyo)<SMALL))
						{
							// Both rods hang out on the same one side
							err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
							if (err==1)
							{
								// New hinge
								hno++;
								(*h).ser=hno;
								(*h).x0=*x;
								(*h).y0=*y;
								(*h).x=0.0L;
								(*h).y=0.0L;
								(*h).rod1=p;
								(*h).rod2=q;
								(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
								h=(*h).next;
								if (!(h==(*hingestart).next)) g=(*g).next;
							}
						}
						else
						{
							if ((fabsl(outx+outxo)<SMALL)&&(fabsl(outy+outyo)<SMALL))
							{
								// The rods hang out only on opposite sides
								err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
								if (err==1)
								{
									// New hinge
									hno++;
									(*h).ser=hno;
									(*h).x0=*x;
									(*h).y0=*y;
									(*h).x=0.0L;
									(*h).y=0.0L;
									(*h).rod1=p;
									(*h).rod2=q;
									(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
									h=(*h).next;
									if (!(h==(*hingestart).next)) g=(*g).next;
								}
							}
							else
							{
								if (((fabsl(outx)<SMALL)&&(fabsl(outyo)<SMALL)&&(fabsl(outxo)>SMALL))||((fabsl(outy)<SMALL)&&(fabsl(outxo)<SMALL)&&(fabsl(outyo)>SMALL)))
								{
									// Second rod hangs out beside first
									err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
									if (err==1)
									{
										// New hinge
										hno++;
										(*h).ser=hno;
										(*h).x0=*x;
										(*h).y0=*y;
										(*h).x=0.0L;
										(*h).y=0.0L;
										(*h).rod1=p;
										(*h).rod2=q;
										(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
										h=(*h).next;
										if (!(h==(*hingestart).next)) g=(*g).next;
									}
									else
									{
										err=intersect((*p).x,(*p).y,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
										if (err==1)
										{
											// New hinge
											hno++;
											(*h).ser=hno;
											(*h).x0=*x;
											(*h).y0=*y;
											(*h).x=0.0L;
											(*h).y=0.0L;
											(*h).rod1=p;
											(*h).rod2=q;
											(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
											h=(*h).next;
											if (!(h==(*hingestart).next)) g=(*g).next;
										}
										else
										{
											err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
											if (err==1)
											{
												// New hinge
												hno++;
												(*h).ser=hno;
												(*h).x0=*x;
												(*h).y0=*y;
												(*h).x=0.0L;
												(*h).y=0.0L;
												(*h).rod1=p;
												(*h).rod2=q;
												(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
												h=(*h).next;
												if (!(h==(*hingestart).next)) g=(*g).next;
											}
											else
											{
												err=intersect((*p).x-outx,(*p).y-outy,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
												if (err==1)
												{
													// New hinge
													hno++;
													(*h).ser=hno;
													(*h).x0=*x;
													(*h).y0=*y;
													(*h).x=0.0L;
													(*h).y=0.0L;
													(*h).rod1=p;
													(*h).rod2=q;
													(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
													h=(*h).next;
													if (!(h==(*hingestart).next)) g=(*g).next;
												}
											}
										}
									}
								}
								else
								{
									// Second rod hangs out on two sides
									if (((fabsl(outx)>SMALL)&&(fabsl(outxo-outx)<SMALL))||((fabsl(outy)>SMALL)&&(fabsl(outyo-outy)<SMALL)))
									{
										// Second rod hangs out the same side as
										// the first, and on another side
										err=intersect((*p).x,(*p).y,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
										if (err==1)
										{
											// New hinge
											hno++;
											(*h).ser=hno;
											(*h).x0=*x;
											(*h).y0=*y;
											(*h).x=0.0L;
											(*h).y=0.0L;
											(*h).rod1=p;
											(*h).rod2=q;
											(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
											h=(*h).next;
											if (!(h==(*hingestart).next)) g=(*g).next;
										}
										else
										{
											err=intersect((*p).x-outx,(*p).y-outy,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
											if (err==1)
											{
												// New hinge
												hno++;
												(*h).ser=hno;
												(*h).x0=*x;
												(*h).y0=*y;
												(*h).x=0.0L;
												(*h).y=0.0L;
												(*h).rod1=p;
												(*h).rod2=q;
												(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
												h=(*h).next;
												if (!(h==(*hingestart).next)) g=(*g).next;
											}
										}
									}
									else
									{
										if (((fabsl(outx)>SMALL)&&(fabsl(outxo+outx)<SMALL))||((fabsl(outy)>SMALL)&&(fabsl(outyo+outy)<SMALL)))
										{
											// Second rod hangs out on the
											// opposite side as the first,
											// and on another side
											err=intersect((*p).x-outx,(*p).y-outy,(*q).x,(*q).y,(*p).tt,(*q).tt,(*p).x0-outx,(*p).y0-outy,(*p).x1-outx,(*p).y1-outy,(*q).x0,(*q).y0,(*q).x1,(*q).y1,x,y);
											if (err==1)
											{
												// New hinge
												hno++;
												(*h).ser=hno;
												(*h).x0=*x;
												(*h).y0=*y;
												(*h).x=0.0L;
												(*h).y=0.0L;
												(*h).rod1=p;
												(*h).rod2=q;
												(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
												h=(*h).next;
												if (!(h==(*hingestart).next)) g=(*g).next;
											}
											else
											{
												err=intersect((*p).x,(*p).y,(*q).x-outxo,(*q).y-outyo,(*p).tt,(*q).tt,(*p).x0,(*p).y0,(*p).x1,(*p).y1,(*q).x0-outxo,(*q).y0-outyo,(*q).x1-outxo,(*q).y1-outyo,x,y);
												if (err==1)
												{
													// New hinge
													hno++;
													(*h).ser=hno;
													(*h).x0=*x;
													(*h).y0=*y;
													(*h).x=0.0L;
													(*h).y=0.0L;
													(*h).rod1=p;
													(*h).rod2=q;
													(*h).next=(struct hinges *)calloc(1,sizeof(struct hinges));
													h=(*h).next;
													if (!(h==(*hingestart).next)) g=(*g).next;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			q=(*q).next;
		}
		p=(*p).next;
	}
	(*g).next=0;
	free(h);
	if (hno==0)
	{
		fprintf(data,"No intersection between rods!!!\n");
		printf("No intersection between rods!!!\n");
		hingestart=0;
	}
	free(x);
	free(y);
	return hno;  // Number of hinges
}

// **************************

int findsticks()
// Find all sticks formed by hinges. Output is number of midpoints (and hence sticks)
{
	struct chain // A chain containing the hinges along one single rod
	{
		struct hinges *h;
		struct chain *next;
	};
	struct rods *r;
	struct hinges *h;
	struct chain *cs,*c,*d,*e;
	struct sticks *s,*t;
	struct midpoints *m,*n;
	int mpno,no,conn;
	long double xd,yd,length,temp;
	mpno=0;
	stickstart=(struct sticks *)calloc(1,sizeof(struct sticks));
	s=stickstart;
	t=s;
	midpointstart=(struct midpoints *)calloc(1,sizeof(struct midpoints));
	m=midpointstart;
	n=m;
	r=rodstart;
	// Go along all rods, and find all hinges along that rod.
	while (r)
	{
		h=hingestart;
		cs=(struct chain *)calloc(1,sizeof(struct chain));
		c=cs;
		d=cs;
//		(*cs).h=0;
		while (h)
		{
			// See if this hinge h is on rod r
			if (((*h).rod1==r)||((*h).rod2==r))
			{
				// Yes, hinge h is on rod r: take note of this hinge
				(*c).h=h;
				(*c).next=(struct chain *)calloc(1,sizeof(struct chain));
				c=(*c).next;
				if (!(c==(*cs).next)) d=(*d).next;
			}
			h=(*h).next;
		}
		// Tie down end of chain: if no hinge along the rod, the chain will already be cleaned up
		(*d).next=0;
		free(c);
		// Now we have all hinges on rod r in chain "chain". Let us order them.
		// Are there hinges along this rod?
		if (!(cs==c))
		{
			// Yes, there is at least one hinge along the rod
			// Is there only one hinge along this rod?
			if (!(cs==d))
			{
				// At least two hinges. First arrange them in order of increasing x0.
				d=cs; // The smallest is to be inserted here, everything is order before this element.
				while (d)
				{
					c=d; // Runs along the chain from d to end to find smallest
					e=d; // This is the smallest we found this far.
					while (c)
					{
						if ((*(*c).h).x0<(*(*e).h).x0)
						{
							e=c;
						}
						c=(*c).next;
					}
					// Now "e" points to the smallest element. Is it the same as d?
					if (e==d)
					{
						// Smallest from d to end is at d: nothing to do
						d=(*d).next;
					}
					else
					{
						// Smallest from d to end is not at d, but later, so it has to be inserted before d
						// Get "e" out of the chain
						c=d;
						while (!((*c).next==e))
						{
							c=(*c).next;
						}
						(*c).next=(*e).next;
						// Now, move "e" before "d"
						// Is "d" the first?
						if (d==cs)
						{
							// d is first, e is new cs
							cs=e;
							(*cs).next=d;
						}
						else
						{
							// d is not first, find element before it
							c=cs;
							while (!((*c).next==d))
							{
								c=(*c).next;
							}
							(*e).next=d;
							(*c).next=e;
						}
					}
				}
				// Now the chain cs is ordered according to x0 (smallest x first).
				// Let us construct the sticks
				c=cs;
				d=(*c).next;
				no=0;   // Number of neighbouring hinge distances larger than 0.5
				conn=0; // Note if first and last hinge form a stick outside domain
				while (d)
				{
					// c and d form a stick inside if they are closer than 0.5 
					// (remember, length of rods <=0.5).
					// Compute distance between c and d:
					xd=fabsl((*(*d).h).x0-(*(*c).h).x0);
					yd=(*(*d).h).y0-(*(*c).h).y0;
					length=sqrtl(xd*xd+yd*yd);
					if (length<0.5)
					{
						// OK, hinges c and d are close, stick is inside
						// Fill in data for stick
						(*t).n1=(*c).h;
						(*t).n2=(*d).h;
						(*t).rod=r;
						(*t).mp=m;
						(*t).len=length;
						// Fill in data for midpoint of stick
						mpno++;
						(*m).ser=mpno;
						(*m).x0=0.5L*((*(*c).h).x0+(*(*d).h).x0);
						(*m).x=0.0L;
						(*m).y0=0.5L*((*(*c).h).y0+(*(*d).h).y0);
						(*m).y=0.0L;
						(*m).n1x=0;
						(*m).n1y=0;
						(*m).n2x=0;
						(*m).n2y=0;
						(*m).st=t;
						(*m).n1=(*c).h;
						(*m).n2=(*d).h;
						// Hinge data on new stick, midpoint and neighbouring hinges
						if (!((*(*c).h).st1))
						{
							(*(*c).h).st1=t;
						}
						else
						{
							if (!((*(*c).h).st2))
							{
								(*(*c).h).st2=t;
							}
							else
							{
								if (!((*(*c).h).st3))
								{
									(*(*c).h).st3=t;
								}
								else
								{
									if (!((*(*c).h).st4))
									{
										(*(*c).h).st4=t;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (c)\n");
									}
								}
							}
						}
						if (!((*(*d).h).st1))
						{
							(*(*d).h).st1=t;
						}
						else
						{
							if (!((*(*d).h).st2))
							{
								(*(*d).h).st2=t;
							}
							else
							{
								if (!((*(*d).h).st3))
								{
									(*(*d).h).st3=t;
								}
								else
								{
									if (!((*(*d).h).st4))
									{
										(*(*d).h).st4=t;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (d)\n");
									}
								}
							}
						}
						if (!((*(*c).h).n1))
						{
							(*(*c).h).n1=(*d).h;
						}
						else
						{
							if (!((*(*c).h).n2))
							{
								(*(*c).h).n2=(*d).h;
							}
							else
							{
								if (!((*(*c).h).n3))
								{
									(*(*c).h).n3=(*d).h;
								}
								else
								{
									if (!((*(*c).h).n4))
									{
										(*(*c).h).n4=(*d).h;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (c)\n");
									}
								}
							}
						}
						if (!((*(*d).h).n1))
						{
							(*(*d).h).n1=(*c).h;
						}
						else
						{
							if (!((*(*d).h).n2))
							{
								(*(*d).h).n2=(*c).h;
							}
							else
							{
								if (!((*(*d).h).n3))
								{
									(*(*d).h).n3=(*c).h;
								}
								else
								{
									if (!((*(*d).h).n4))
									{
										(*(*d).h).n4=(*c).h;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (d)\n");
									}
								}
							}
						}
						if (!((*(*c).h).mp1))
						{
							(*(*c).h).mp1=m;
						}
						else
						{
							if (!((*(*c).h).mp2))
							{
								(*(*c).h).mp2=m;
							}
							else
							{
								if (!((*(*c).h).mp3))
								{
									(*(*c).h).mp3=m;
								}
								else
								{
									if (!((*(*c).h).mp4))
									{
										(*(*c).h).mp4=m;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (c)\n");
									}
								}
							}
						}
						if (!((*(*d).h).mp1))
						{
							(*(*d).h).mp1=m;
						}
						else
						{
							if (!((*(*d).h).mp2))
							{
								(*(*d).h).mp2=m;
							}
							else
							{
								if (!((*(*d).h).mp3))
								{
									(*(*d).h).mp3=m;
								}
								else
								{
									if (!((*(*d).h).mp4))
									{
										(*(*d).h).mp4=m;
									}
									else
									{
										printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (d)\n");
									}
								}
							}
						}
						// Room for next stick and midpoint
						(*t).next=(struct sticks *)calloc(1,sizeof(struct sticks));
						s=t;
						t=(*t).next;
						(*m).next=(struct midpoints *)calloc(1,sizeof(struct midpoints));
						n=m;
						m=(*m).next;
					}
					else
					{
						// Now the neighbouring hinges are further than 0.5.
						no++;
						// This implies there might be a stick partially outside:
						// check if it is the case.
						if (no==1)
						{
							// First gap between hinges larger than 0.5
							if (xd<0.5)
							{
								// Hanging out on top and bottom: there is a stick here
								// Fill in data for stick
								(*t).n1=(*c).h;
								(*t).n2=(*d).h;
								(*t).rod=r;
								(*t).mp=m;
								// Compute real length
								if (yd<0.0L)
								{
									length=yd+WY;
								}
								else
								{
									length=WY-yd;
								}
								(*t).len=sqrtl(length*length+xd*xd);
								// Fill in data for midpoint of stick
								mpno++;
								(*m).n1x=0;
								(*m).n1y=0;
								(*m).n2x=0;
								(*m).n2y=0;
								(*m).ser=mpno;
								(*m).x0=0.5L*((*(*c).h).x0+(*(*d).h).x0);
								(*m).x=0.0L;
								temp=0.5L*((*(*c).h).y0+(*(*d).h).y0+WY);
								if (temp<WY)
								{
									if ((*(*c).h).y0 < (*(*d).h).y0)
									{
										(*m).n1y=1;
									}
									else
									{
										(*m).n2y=1;
									}
								}
								else
								{
									temp=temp-WY;
									if ((*(*c).h).y0 < (*(*d).h).y0)
									{
										(*m).n2y=-1;
									}
									else
									{
										(*m).n1y=-1;
									}
		
								}
								(*m).y0=temp;;
								(*m).y=0.0L;
								(*m).st=t;
								(*m).n1=(*c).h;
								(*m).n2=(*d).h;
								// Hinge data on new stick, midpoint and neighbouring hinges
								if (!((*(*c).h).st1))
								{
									(*(*c).h).st1=t;
								}
								else
								{
									if (!((*(*c).h).st2))
									{
										(*(*c).h).st2=t;
									}
									else
									{
										if (!((*(*c).h).st3))
										{
											(*(*c).h).st3=t;
										}
										else
										{
											if (!((*(*c).h).st4))
											{
												(*(*c).h).st4=t;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).st1))
								{
									(*(*d).h).st1=t;
								}
								else
								{
									if (!((*(*d).h).st2))
									{
										(*(*d).h).st2=t;
									}
									else
									{
										if (!((*(*d).h).st3))
										{
											(*(*d).h).st3=t;
										}
										else
										{
											if (!((*(*d).h).st4))
											{
												(*(*d).h).st4=t;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (d)\n");
											}
										}
									}
								}
								if (!((*(*c).h).n1))
								{
									(*(*c).h).n1=(*d).h;
								}
								else
								{
									if (!((*(*c).h).n2))
									{
										(*(*c).h).n2=(*d).h;
									}
									else
									{
										if (!((*(*c).h).n3))
										{
											(*(*c).h).n3=(*d).h;
										}
										else
										{
											if (!((*(*c).h).n4))
											{
												(*(*c).h).n4=(*d).h;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).n1))
								{
									(*(*d).h).n1=(*c).h;
								}
								else
								{
									if (!((*(*d).h).n2))
									{
										(*(*d).h).n2=(*c).h;
									}
									else
									{
										if (!((*(*d).h).n3))
										{
											(*(*d).h).n3=(*c).h;
										}
										else
										{
											if (!((*(*d).h).n4))
											{
												(*(*d).h).n4=(*c).h;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (d)\n");
											}
										}
									}
								}
								if (!((*(*c).h).mp1))
								{
									(*(*c).h).mp1=m;
								}
								else
								{
									if (!((*(*c).h).mp2))
									{
										(*(*c).h).mp2=m;
									}
									else
									{
										if (!((*(*c).h).mp3))
										{
											(*(*c).h).mp3=m;
										}
										else
										{
											if (!((*(*c).h).mp4))
											{
												(*(*c).h).mp4=m;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).mp1))
								{
									(*(*d).h).mp1=m;
								}
								else
								{
									if (!((*(*d).h).mp2))
									{
										(*(*d).h).mp2=m;
									}
									else
									{
										if (!((*(*d).h).mp3))
										{
											(*(*d).h).mp3=m;
										}
										else
										{
											if (!((*(*d).h).mp4))
											{
												(*(*d).h).mp4=m;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (d)\n");
											}
										}
									}
								}
								// Room for next stick and midpoint
								(*t).next=(struct sticks *)calloc(1,sizeof(struct sticks));
								s=t;
								t=(*t).next;
								(*m).next=(struct midpoints *)calloc(1,sizeof(struct midpoints));
								n=m;
								m=(*m).next;
							}
							else
							{
								// Hanging out on left and right: no stick.
								// However, note that first and last hinge are connected by stick
								conn=1;
							}
						}
						if (no==2)
						{
							// Second gap between hinges larger than 0.5
							if (xd<0.5)
							{
								// Hanging out on top and bottom: there is a stick here
								// Fill in data for stick
								(*t).n1=(*c).h;
								(*t).n2=(*d).h;
								(*t).rod=r;
								(*t).mp=m;
								// Compute real length
								if (yd<0.0L)
								{
									length=yd+WY;
								}
								else
								{
									length=WY-yd;
								}
								(*t).len=sqrtl(length*length+xd*xd);
								// Fill in data for midpoint of stick
								mpno++;
								(*m).n1x=0;
								(*m).n1y=0;
								(*m).n2x=0;
								(*m).n2y=0;
								(*m).ser=mpno;
								(*m).x0=0.5L*((*(*c).h).x0+(*(*d).h).x0);
								(*m).x=0.0L;
								temp=0.5L*((*(*c).h).y0+(*(*d).h).y0+WY);
								if (temp<WY)
								{
									if ((*(*c).h).y0 < (*(*d).h).y0)
									{
										(*m).n1y=1;
									}
									else
									{
										(*m).n2y=1;
									}
								}
								else
								{
									temp=temp-WY;
									if ((*(*c).h).y0 < (*(*d).h).y0)
									{
										(*m).n2y=-1;
									}
									else
									{
										(*m).n1y=-1;
									}
		
								}
								(*m).y0=temp;;
								(*m).y=0.0L;
								(*m).st=t;
								(*m).n1=(*c).h;
								(*m).n2=(*d).h;
								// Hinge data on new stick, midpoint and neighbouring hinges
								if (!((*(*c).h).st1))
								{
									(*(*c).h).st1=t;
								}
								else
								{
									if (!((*(*c).h).st2))
									{
										(*(*c).h).st2=t;
									}
									else
									{
										if (!((*(*c).h).st3))
										{
											(*(*c).h).st3=t;
										}
										else
										{
											if (!((*(*c).h).st4))
											{
												(*(*c).h).st4=t;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).st1))
								{
									(*(*d).h).st1=t;
								}
								else
								{
									if (!((*(*d).h).st2))
									{
										(*(*d).h).st2=t;
									}
									else
									{
										if (!((*(*d).h).st3))
										{
											(*(*d).h).st3=t;
										}
										else
										{
											if (!((*(*d).h).st4))
											{
												(*(*d).h).st4=t;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (d)\n");
											}
										}
									}
								}
								if (!((*(*c).h).n1))
								{
									(*(*c).h).n1=(*d).h;
								}
								else
								{
									if (!((*(*c).h).n2))
									{
										(*(*c).h).n2=(*d).h;
									}
									else
									{
										if (!((*(*c).h).n3))
										{
											(*(*c).h).n3=(*d).h;
										}
										else
										{
											if (!((*(*c).h).n4))
											{
												(*(*c).h).n4=(*d).h;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).n1))
								{
									(*(*d).h).n1=(*c).h;
								}
								else
								{
									if (!((*(*d).h).n2))
									{
										(*(*d).h).n2=(*c).h;
									}
									else
									{
										if (!((*(*d).h).n3))
										{
											(*(*d).h).n3=(*c).h;
										}
										else
										{
											if (!((*(*d).h).n4))
											{
												(*(*d).h).n4=(*c).h;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (d)\n");
											}
										}
									}
								}
								if (!((*(*c).h).mp1))
								{
									(*(*c).h).mp1=m;
								}
								else
								{
									if (!((*(*c).h).mp2))
									{
										(*(*c).h).mp2=m;
									}
									else
									{
										if (!((*(*c).h).mp3))
										{
											(*(*c).h).mp3=m;
										}
										else
										{
											if (!((*(*c).h).mp4))
											{
												(*(*c).h).mp4=m;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (c)\n");
											}
										}
									}
								}
								if (!((*(*d).h).mp1))
								{
									(*(*d).h).mp1=m;
								}
								else
								{
									if (!((*(*d).h).mp2))
									{
										(*(*d).h).mp2=m;
									}
									else
									{
										if (!((*(*d).h).mp3))
										{
											(*(*d).h).mp3=m;
										}
										else
										{
											if (!((*(*d).h).mp4))
											{
												(*(*d).h).mp4=m;
											}
											else
											{
												printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (d)\n");
											}
										}
									}
								}
								// Room for next stick and midpoint
								(*t).next=(struct sticks *)calloc(1,sizeof(struct sticks));
								s=t;
								t=(*t).next;
								(*m).next=(struct midpoints *)calloc(1,sizeof(struct midpoints));
								n=m;
								m=(*m).next;
							}
							else
							{
								// Hanging out on left and right: no stick.
								// However, note that first and last hinge are connected by stick
								conn=1;
							}

						}
						if (no>2)
						{
							printf("ERROR: ROD CANNOT HANG OUT MORE THAN TWICE!!!\n");
						}
					}
					// Next two hinges
					d=(*d).next;
					c=(*c).next;
				}
				// Check if first and last hinges are connected by a stick
				if (conn)
				{
					// Yes, a stick connects first and last hinge outside the domain
					// Rename hinges
					d=c;
					c=cs;
					// Fill in data for stick
					(*t).n1=(*c).h;
					(*t).n2=(*d).h;
					(*t).rod=r;
					(*t).mp=m;
					// Compute real length
					xd=(*(*c).h).x0+WX-(*(*d).h).x0;
					yd=(*(*c).h).y0-(*(*d).h).y0;
					if (fabsl(yd)>0.5L)
					{
						if (yd<0.0L)
						{
							yd=yd+WY;
						}
						else
						{
							yd=WY-yd;
						}
					}
					(*t).len=sqrtl(xd*xd+yd*yd);
					// Fill in data for midpoint of stick
					mpno++;
					(*m).n1x=0;
					(*m).n1y=0;
					(*m).n2x=0;
					(*m).n2y=0;
					(*m).ser=mpno;
					temp=0.5L*((*(*c).h).x0+(*(*d).h).x0+WX);
					if (temp<WX)
					{
						(*m).n1x=1;
					}
					else
					{
						temp=temp-WX;
						(*m).n2x=-1;
					}
					(*m).x0=temp;
					(*m).x=0.0L;
					yd=(*(*c).h).y0-(*(*d).h).y0;
					if (fabsl(yd)>0.5L)
					{
						temp=0.5L*((*(*c).h).y0+(*(*d).h).y0+WY);
						if (temp<WY)
						{
							if ((*(*c).h).y0<(*(*d).h).y0)
							{
								(*m).n1y=1;
							}
							else
							{
								(*m).n2y=1;
							}
						}
						else
						{
							temp=temp-WY;
							if ((*(*c).h).y0<(*(*d).h).y0)
							{
								(*m).n2y=-1;
							}
							else
							{
								(*m).n1y=-1;
							}
						}
						(*m).y0=temp;
					}
					else
					{
						(*m).y0=0.5L*((*(*c).h).y0+(*(*d).h).y0);
					}
					(*m).y=0.0L;
					(*m).st=t;
					(*m).n1=(*c).h;
					(*m).n2=(*d).h;
					// Hinge data on new stick, midpoint and neighbouring hinges
					if (!((*(*c).h).st1))
					{
						(*(*c).h).st1=t;
					}
					else
					{
						if (!((*(*c).h).st2))
						{
							(*(*c).h).st2=t;
						}
						else
						{
							if (!((*(*c).h).st3))
							{
								(*(*c).h).st3=t;
							}
							else
							{
								if (!((*(*c).h).st4))
								{
									(*(*c).h).st4=t;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (c)\n");
								}
							}
						}
					}
					if (!((*(*d).h).st1))
					{
						(*(*d).h).st1=t;
					}
					else
					{
						if (!((*(*d).h).st2))
						{
							(*(*d).h).st2=t;
						}
						else
						{
							if (!((*(*d).h).st3))
							{
								(*(*d).h).st3=t;
							}
							else
							{
								if (!((*(*d).h).st4))
								{
									(*(*d).h).st4=t;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 STICKS (d)\n");
								}
							}
						}
					}
					if (!((*(*c).h).n1))
					{
						(*(*c).h).n1=(*d).h;
					}
					else
					{
						if (!((*(*c).h).n2))
						{
							(*(*c).h).n2=(*d).h;
						}
						else
						{
							if (!((*(*c).h).n3))
							{
								(*(*c).h).n3=(*d).h;
							}
							else
							{
								if (!((*(*c).h).n4))
								{
									(*(*c).h).n4=(*d).h;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (c)\n");
								}
							}
						}
					}
					if (!((*(*d).h).n1))
					{
						(*(*d).h).n1=(*c).h;
					}
					else
					{
						if (!((*(*d).h).n2))
						{
							(*(*d).h).n2=(*c).h;
						}
						else
						{
							if (!((*(*d).h).n3))
							{
								(*(*d).h).n3=(*c).h;
							}
							else
							{
								if (!((*(*d).h).n4))
								{
									(*(*d).h).n4=(*c).h;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 NEIGHBOURS (d)\n");
								}
							}
						}
					}
					if (!((*(*c).h).mp1))
					{
						(*(*c).h).mp1=m;
					}
					else
					{
						if (!((*(*c).h).mp2))
						{
							(*(*c).h).mp2=m;
						}
						else
						{
							if (!((*(*c).h).mp3))
							{
								(*(*c).h).mp3=m;
							}
							else
							{
								if (!((*(*c).h).mp4))
								{
									(*(*c).h).mp4=m;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (c)\n");
								}
							}
						}
					}
					if (!((*(*d).h).mp1))
					{
						(*(*d).h).mp1=m;
					}
					else
					{
						if (!((*(*d).h).mp2))
						{
							(*(*d).h).mp2=m;
						}
						else
						{
							if (!((*(*d).h).mp3))
							{
								(*(*d).h).mp3=m;
							}
							else
							{
								if (!((*(*d).h).mp4))
								{
									(*(*d).h).mp4=m;
								}
								else
								{
									printf("ERROR: NOT POSSIBLE TO HAVE 5 MIDPOINTS (d)\n");
								}
							}
						}
					}
					// Room for next stick and midpoint
					(*t).next=(struct sticks *)calloc(1,sizeof(struct sticks));
					s=t;
					t=(*t).next;
					(*m).next=(struct midpoints *)calloc(1,sizeof(struct midpoints));
					n=m;
					m=(*m).next;
				}
			}
			// Clean up "chain"
			c=cs;
			while (c)
			{
				c=(*cs).next;
				free(cs);
				cs=c;
			}
		}
		// This rod is done, all sticks are found, next rod
		r=(*r).next;
	}
	// Clean last element of stick, midpoint
	(*s).next=0;
	free(t);
	(*n).next=0;
	free(m);
	return mpno;
}

// **************************

int deletehinges()
// Deletes hinges with no attached sticks.
{
	int hno;
	struct hinges *h,*hn;
	hno=0;
	// See if first hinge(s) has(ve) attached stick
	while (!((*hingestart).st1))
	{
		h=hingestart;
		hingestart=(*hingestart).next;
		free(h);
		if (!(hingestart))
		{
			// If this was the last hinge, jump out
			break;
		}
	}
	// Do we have hinges remaining?
	if (hingestart)
	{
		// Set serial number of stick
		hno++;
		(*hingestart).ser=1;
		h=hingestart;
		hn=(*hingestart).next;
		while (hn)
		{
			if ((*hn).st1)
			{
				hno++;
				(*hn).ser=hno;
				h=hn;
				hn=(*hn).next;
			}
			else
			{
				// Actual hinge is dangling
				(*h).next=(*hn).next;
				free(hn);
				hn=(*h).next;
			}
		}
	}
	return hno;
}

// **************************

int printrods()
// Prints all rods
{
	struct rods *p;
	p=rodstart;
	while (p)
	{
		fprintf(roddata,"%20p  %Lf  %Lf  %9.6Lf  %12.6Lf  %9.6Lf  %9.6Lf  %9.6Lf  %9.6Lf  %9.6Lf  %9.6Lf\n",(void *)p,(*p).x,(*p).y,(*p).teta,(*p).tt,(*p).sit,(*p).cot,(*p).x0,(*p).y0,(*p).x1,(*p).y1);
		p=(*p).next;
	}
	fprintf(roddata,"\n\n");
	return 0;
}

// **************************

int printhinges()
// Prints all hinges
{
	struct hinges *h;
	h=hingestart;
	while (h)
	{
		fprintf(hingedata,"%20p  %10ld  %Lf  %Lf  %Lf  %Lf  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p  %20p\n",(void *)h,(*h).ser,(*h).x0,(*h).y0,(*h).x,(*h).y,(void *)(*h).rod1,(void *)(*h).rod2,(void *)(*h).st1,(void *)(*h).st2,(void *)(*h).st3,(void *)(*h).st4,(void *)(*h).n1,(void *)(*h).n2,(void *)(*h).n3,(void *)(*h).n4,(void *)(*h).mp1,(void *)(*h).mp2,(void *)(*h).mp3,(void *)(*h).mp4);
		h=(*h).next;
	}
	fprintf(hingedata,"\n\n");
	return 0;
}

// **************************

int printmidpoints()
// Prints all midpoints
{
	struct midpoints *h;
	h=midpointstart;
	while (h)
	{
		fprintf(midpointdata,"%20p  %10ld  %Lf  %Lf  %Lf  %Lf  %20p  %20p  %20p  %3d  %3d  %3d  %3d\n",(void *)h,(*h).ser,(*h).x0,(*h).y0,(*h).x,(*h).y,(void *)(*h).st,(void *)(*h).n1,(void *)(*h).n2,(*h).n1x,(*h).n1y,(*h).n2x,(*h).n2y);
		h=(*h).next;
	}
	fprintf(midpointdata,"\n\n");
	return 0;
}

// **************************

int printsticks(long double shl, long double stl)
// Prints all stick data
// Input: shearload (shl) and stretchload (stl),
// both are displacements.
{
	struct sticks *h;
	long double x1,y1,x2,y2;
	h=stickstart;
	while (h)
	{
		x1=0;
		y1=0;
		x2=0;
		y2=0;
		if ((*(*h).mp).n1x==-1)
		{
			x1=x1-WX-stl;
			y1=y1-shl;
		}
		if ((*(*h).mp).n1x==1)
		{
			x1=x1+WX+stl;
			y1=y1+shl;
		}
		if ((*(*h).mp).n1y==-1)
		{
			y1=y1-WY;
		}
		if ((*(*h).mp).n1y==1)
		{
			y1=y1+WY;
		}
		if ((*(*h).mp).n2x==-1)
		{
			x2=x2-WX-stl;
			y2=y2-shl;
		}
		if ((*(*h).mp).n2x==1)
		{
			x2=x2+WX+stl;
			y2=y2+shl;
		}
		if ((*(*h).mp).n2y==-1)
		{
			y2=y2-WY;
		}
		if ((*(*h).mp).n2y==1)
		{
			y2=y2+WY;
		}
		x1=x1+(*(*h).n1).x+(*(*h).n1).x0;
		y1=y1+(*(*h).n1).y+(*(*h).n1).y0;
		x2=x2+(*(*h).n2).x+(*(*h).n2).x0;
		y2=y2+(*(*h).n2).y+(*(*h).n2).y0;
		fprintf(stickdata,"%20p  %20p  %20p  %20p  %20p  %Lf  %Lf  %Lf  %Lf  %Lf\n",(void *)h,(void *)(*h).n1,(void *)(*h).n2,(void *)(*h).mp,(void *)(*h).rod,(*h).len,x1,y1,x2,y2);
		h=(*h).next;
	}
	fprintf(stickdata,"\n\n");
	return 0;
}

// **************************

int setupfiles()
// Opens files for writing and prints header lines into files, and also basic data
{
	char fn[30];
	// Open output files
	sprintf(fn,"%s%04d.dat",DATAFILE,FILENUMBER);
	data=fopen(fn,"w");
	sprintf(fn,"%s%04d.dat",RODFILE,FILENUMBER);
	roddata=fopen(fn,"w");
	sprintf(fn,"%s%04d.dat",HINGEFILE,FILENUMBER);
	hingedata=fopen(fn,"w");
	sprintf(fn,"%s%04d.dat",MIDPOINTFILE,FILENUMBER);
	midpointdata=fopen(fn,"w");
	sprintf(fn,"%s%04d.dat",STICKFILE,FILENUMBER);
	stickdata=fopen(fn,"w");
	// Print some parameters
	fprintf(data,"# WX=%Lf  WY=%Lf  L=%Lf  N=%d  SEED=%d  CASENO=%d\n",WX,WY,L,N,SEED,CASENO);
	fprintf(data,"# EA=%.12Lf  EI=%.12Lf  GI=%.12Lf\n",EA,EI,GI);
	fprintf(data,"# SHEARDISPL=%.15Lf  STRETCHDISPL=%.15Lf\n",SHEARDISPL,STRETCHDISPL);
	// Print columns for datafiles
	fprintf(roddata,"#       p->rod          mid x     mid y     angle        tangent      sin        cos      start x    start y     end x      end y\n");
	fprintf(hingedata,"#       p->hinge         ser.no.   orig x    orig y    delta x   delta y            p->rod1               p->rod2               p->st1                p->st2                p->st3                p->st4                p->nh1                p->nh2                p->nh3                p->nh4                p->mp1                p->mp2                p->mp3                p->mp4\n");
	fprintf(midpointdata,"#         p->mp          ser.no.   orig x    orig y    delta x   delta y            p->st                 p->nh1                p->nh2       n1x  n1y  n2x  n2y\n");
	fprintf(stickdata,"#       p->st                   p->n1                 p->n2                 p->mp                 p->rod       length      x0        y0        x1        y1\n");
	fprintf(data,"#        lc                 lb                 lambda             L/Lc               L/lambda           lb/L          E_stretch            E_bend               E_shear              E_total              E_fromA              (E_total-E_fromA)/E_total\n");
	return 0;
}


// **************************

int networkdata()
// Prints the basic data of the actual network
{
	struct sticks *s;
	int err;
	long double lc,lb,lambda,lplc,lpla,lbpl;
	// Stiffness ratio (m)
	lb=sqrtl(EI/EA);
	// Average distance between hinges (lc) (m)
	err=0;
	lc=0.0L;
	s=stickstart;
	while (s)
	{
		err++;
		lc=lc+(*s).len;
		s=(*s).next;
	}
	lc=lc/(long double)err;
	// Another distance (m)
	lambda=powl(lc*lc*lc*lc/lb,1.0L/3.0L);
	// Dimensionless numbers
	lplc=L/lc;  // Density of net: mean number of hinges on a rod
	lpla=L/lambda;
	lbpl=lb/L;
	// Print parameters
	fprintf(data,"%18.9Lf %18.9Lf %18.9Lf %18.9Lf %18.9Lf %18.9Lf ", lc, lb, lambda, lplc, lpla, lbpl);
	return 0;
}

// **************************

int cleanuprods()
// Cleans up all rods
{
	struct rods *p;
	p=rodstart;
	while (p)
	{
		p=(*rodstart).next;
		free(rodstart);
		rodstart=p;
	}
	return 0;
}

// **************************

int cleanuphinges()
// Cleans up all hinges
{
	struct hinges *p;
	p=hingestart;
	while (p)
	{
		p=(*hingestart).next;
		free(hingestart);
		hingestart=p;
	}
	return 0;
}

// **************************

int cleanupmidpoints()
// Cleans up all midpoints
{
	struct midpoints *p;
	p=midpointstart;
	while (p)
	{
		p=(*midpointstart).next;
		free(midpointstart);
		midpointstart=p;
	}
	return 0;
}

// **************************

int cleanupsticks()
// Cleans up all sticks
{
	struct sticks *p;
	p=stickstart;
	while (p)
	{
		p=(*stickstart).next;
		free(stickstart);
		stickstart=p;
	}
	return 0;
}

// **************************

long double stretchenergy(long double shl, long double stl)
// Computes the energy stored in the stretch of all half sticks
{
	long double energy,x0,y0,xa,ya,l0,dum,st,ct;
	struct midpoints *p;
	// Cycle along all midpoints
	energy=0.0;
	p=midpointstart;
	while (p)
	{
		// Length of half stick, sine and cosine of angle
		l0=0.5L*(*(*p).st).len;
		st=(*(*(*p).st).rod).sit;
		ct=(*(*(*p).st).rod).cot;
		// Nodal displacements: midpoint
		x0=(*p).x;
		y0=(*p).y;
		// Nodal displacements: first hinge, 
		// modified by load displacement if hinge hangs out
		xa=(*(*p).n1).x+stl*((long double)(*p).n1x);
		ya=(*(*p).n1).y+shl*((long double)(*p).n1x);
		dum=ct*(x0-xa)+st*(y0-ya);
		energy=energy+dum*dum*EA/(2.0L*l0);
		// Nodal displacements: second hinge, 
		// modified by load displacement if hinge hangs out
		xa=(*(*p).n2).x+stl*((long double)(*p).n2x);
		ya=(*(*p).n2).y+shl*((long double)(*p).n2x);
		dum=ct*(x0-xa)+st*(y0-ya);
		energy=energy+dum*dum*EA/(2.0L*l0);
		p=(*p).next;
	}
	return energy;
}

// **************************

long double bendenergy(long double shl, long double stl)
// Computes the energy stored in the bending of all sticks
{
	int cont1,cont2,n1x;
	long double energy,sg,cg,l1,l2,dum,xi,yi,xj,yj,xk,yk;
	struct midpoints *p;
	struct hinges *h;
	// Cycle along all midpoints to compute bend of center of sticks
	energy=0.0;
	p=midpointstart;
	while (p)
	{
		// Compute bend of actual stick
		l1=0.5L*(*(*p).st).len;
		sg=(*(*(*p).st).rod).sit;  // Sine of rod angle
		cg=(*(*(*p).st).rod).cot;  // Cosine of rod angle
		// Nodal displacements: midpoint
		xj=(*p).x;
		yj=(*p).y;
		// Nodal displacements: first hinge, 
		// modified by load displacement if hinge hangs out
		xi=(*(*p).n1).x+stl*((long double)(*p).n1x);
		yi=(*(*p).n1).y+shl*((long double)(*p).n1x);
		// Nodal displacements: second hinge, 
		// modified by load displacement if hinge hangs out
		xk=(*(*p).n2).x+stl*((long double)(*p).n2x);
		yk=(*(*p).n2).y+shl*((long double)(*p).n2x);
		// Contribution to energy
		dum=sg*(xi+xk-2.0L*xj)-cg*(yi+yk-2.0L*yj);
		dum=EI*dum*dum/(2.0L*l1*l1*l1);
		energy=energy+dum;
		// Check if there is a stick on the other side 
		// of the neighbouring hinges: if not then energy is increased
		// so that later we need to consider only hinges with sticks
		// on both sides.
		// Hinge on one end:
		cont1=1;
		if (((*(*p).n1).st1)==((*p).st))
		{
			if (((*(*p).n1).st2)==0)
			{
				cont1=0;
			}
			else
			{
				if (((*(*(*p).n1).st2).rod)!=((*(*p).st).rod))
				{
					cont1=0;
				}
			}
		}
		if (((*(*p).n1).st2)==((*p).st))
		{
			if (((*(*(*p).n1).st1).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n1).st3)==0)
				{
					cont1=0;
				}
				else
				{
					if (((*(*(*p).n1).st3).rod)!=((*(*p).st).rod))
					{
						printf("(A0) ILYEN ESET NEM LEHET!!!\n");
						cont1=0;
					}
				}
			}
		}
		if (((*(*p).n1).st3)==((*p).st))
		{
			if (((*(*(*p).n1).st2).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n1).st4)==0)
				{
					cont1=0;
				}
				else
				{
					if (((*(*(*p).n1).st4).rod)!=((*(*p).st).rod))
					{
						printf("(B0) ILYEN ESET NEM LEHET!!!\n");
						cont1=0;
					}
				}
			}
		}
		if (((*(*p).n1).st4)==((*p).st))
		{
			if (((*(*(*p).n1).st3).rod)!=((*(*p).st).rod))
			{
				cont1=0;
			}
		}
		if (cont1==0)
		{
			energy=energy+0.5L*dum;
		}
		// Hinge on other end:
		cont2=1;
		if (((*(*p).n2).st1)==((*p).st))
		{
			if (((*(*p).n2).st2)==0)
			{
				cont2=0;
			}
			else
			{
				if (((*(*(*p).n2).st2).rod)!=((*(*p).st).rod))
				{
					cont2=0;
				}
			}
		}
		if (((*(*p).n2).st2)==((*p).st))
		{
			if (((*(*(*p).n2).st1).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n2).st3)==0)
				{
					cont2=0;
				}
				else
				{
					if (((*(*(*p).n2).st3).rod)!=((*(*p).st).rod))
					{
						printf("(C0) ILYEN ESET NEM LEHET!!!\n");
						cont2=0;
					}
				}
			}
		}
		if (((*(*p).n2).st3)==((*p).st))
		{
			if (((*(*(*p).n2).st2).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n2).st4)==0)
				{
					cont2=0;
				}
				else
				{
					if (((*(*(*p).n2).st4).rod)!=((*(*p).st).rod))
					{
						printf("(D0) ILYEN ESET NEM LEHET!!!\n");
						cont2=0;
					}
				}
			}
		}
		if (((*(*p).n2).st4)==((*p).st))
		{
			if (((*(*(*p).n2).st3).rod)!=((*(*p).st).rod))
			{
				cont2=0;
			}
		}
		if (cont2==0)
		{
			energy=energy+0.5L*dum;
		}
		// Next midpoint
		p=(*p).next;
	}
	// Cycle along all hinges to compute bend of joining sticks
	h=hingestart;
	while (h)
	{
		// Compute bend of sticks at actual hinge
		// Need only to consider sticks if there is a stick on
		// both sides of the hinge, dangling ends are treated already

		// Check if there is any stick at all
		if (((*h).st1)&&((*h).st2))
		{
			// There are two sticks, check if they belong to same rod
			if (((*(*h).st1).rod)==((*(*h).st2).rod))
			{
				// OK, sticks 1 and 2 are on the same rod, compute energy
				l1=0.5L*(*(*h).st1).len;
				l2=0.5L*(*(*h).st2).len;
				sg=(*(*(*(*h).mp1).st).rod).sit;
				cg=(*(*(*(*h).mp1).st).rod).cot;
				// Nodal displacements: hinge
				xj=(*h).x;
				yj=(*h).y;
				// Nodal displacements: first midpoint,
				xi=(*(*h).mp1).x;
				yi=(*(*h).mp1).y;
				// Nodal displacements: second midpoint,
				xk=(*(*h).mp2).x;
				yk=(*(*h).mp2).y;
				// Midpoint displacements are modified by load if hinge hangs out
				// First midpoint
				if ((*(*h).mp1).n1==h)
				{
					n1x=(*(*h).mp1).n1x;
				}
				else
				{
					n1x=(*(*h).mp1).n2x;
				}
				xi=xi-((long double)n1x)*stl;
				yi=yi-((long double)n1x)*shl;
				// Second midpoint
				if ((*(*h).mp2).n1==h)
				{
					n1x=(*(*h).mp2).n1x;
				}
				else
				{
					n1x=(*(*h).mp2).n2x;
				}
				xk=xk-((long double)n1x)*stl;
				yk=yk-((long double)n1x)*shl;
				dum=sg*((xi-xj)/l1+(xk-xj)/l2)+cg*((yj-yi)/l1+(yj-yk)/l2);
				energy=energy+EI*dum*dum/(l1+l2);
				// Check if there are more sticks
				if (((*h).st3)&&((*h).st4))
				{
					// There are two more sticks, are they on the same rod?
					if (((*(*h).st3).rod)==((*(*h).st4).rod))
					{
						// OK, sticks 3 and 4 are on the same rod, compute energy
						l1=0.5L*(*(*h).st3).len;
						l2=0.5L*(*(*h).st4).len;
						sg=(*(*(*(*h).mp3).st).rod).sit;
						cg=(*(*(*(*h).mp3).st).rod).cot;
						// Nodal displacements: hinge
						xj=(*h).x;
						yj=(*h).y;
						// Nodal displacements: first midpoint,
						xi=(*(*h).mp3).x;
						yi=(*(*h).mp3).y;
						// Nodal displacements: second midpoint,
						xk=(*(*h).mp4).x;
						yk=(*(*h).mp4).y;
						// Midpoint displacements are modified by load if hinge hangs out
						// First midpoint
						if ((*(*h).mp3).n1==h)
						{
							n1x=(*(*h).mp3).n1x;
						}
						else
						{
							n1x=(*(*h).mp3).n2x;
						}
						xi=xi-((long double)n1x)*stl;
						yi=yi-((long double)n1x)*shl;
						// Second midpoint
						if ((*(*h).mp4).n1==h)
						{
							n1x=(*(*h).mp4).n1x;
						}
						else
						{
							n1x=(*(*h).mp4).n2x;
						}
						xk=xk-((long double)n1x)*stl;
						yk=yk-((long double)n1x)*shl;
						dum=sg*((xi-xj)/l1+(xk-xj)/l2)+cg*((yj-yi)/l1+(yj-yk)/l2);
						energy=energy+EI*dum*dum/(l1+l2);
					}
				}
      			}
			else
			{
				// So first two sticks not on the same rod, are there more sticks?
				if ((*h).st3)
				{
					// Are sticks 2 and 3 on the same rod?
					if (((*(*h).st2).rod)==((*(*h).st3).rod))
					{
						// OK, sticks 2 and 3 are on the same rod, compute energy
						l1=0.5L*(*(*h).st2).len;
						l2=0.5L*(*(*h).st3).len;
						sg=(*(*(*(*h).mp2).st).rod).sit;
						cg=(*(*(*(*h).mp2).st).rod).cot;
						// Nodal displacements: hinge
						xj=(*h).x;
						yj=(*h).y;
						// Nodal displacements: first midpoint,
						xi=(*(*h).mp2).x;
						yi=(*(*h).mp2).y;
						// Nodal displacements: second midpoint,
						xk=(*(*h).mp3).x;
						yk=(*(*h).mp3).y;
						// Midpoint displacements are modified by load if hinge hangs out
						// First midpoint
						if ((*(*h).mp2).n1==h)
						{
							n1x=(*(*h).mp2).n1x;
						}
						else
						{
							n1x=(*(*h).mp2).n2x;
						}
						xi=xi-((long double)n1x)*stl;
						yi=yi-((long double)n1x)*shl;
						// Second midpoint
						if ((*(*h).mp3).n1==h)
						{
							n1x=(*(*h).mp3).n1x;
						}
						else
						{
							n1x=(*(*h).mp3).n2x;
						}
						xk=xk-((long double)n1x)*stl;
						yk=yk-((long double)n1x)*shl;
						dum=sg*((xi-xj)/l1+(xk-xj)/l2)+cg*((yj-yi)/l1+(yj-yk)/l2);
						energy=energy+EI*dum*dum/(l1+l2);
					}
				}

			}
		}
		// Next hinge
		h=(*h).next;
	}
	return energy;
}

// **************************

long double shearenergy(long double shl, long double stl)
// Computes the energy stored in the shear at all hinges
{
	struct hinges *h;
	long double energy,nxj,nxk,nxl,nxm,xj,yj,xk,yk,xl,yl,xm,ym,lj,lk,ll,lm,sg1,cg1,sg2,cg2,dg;
	// Cycle along all hinges to compute energy from shear springs
	h=hingestart;
	energy=0.0L;
	while (h)
	{
		// Check how many sticks join here, need to have two connecting rods here with valid sticks
		if (((*h).st3) || (((*h).st2)&&(((*(*h).st1).rod)!=((*(*h).st2).rod))))
		{
			// We have at least 3 sticks or two sticks on different rods
			// Hence we have shear energy at this hinge
			// Common for all cases:
			if ((*(*h).mp1).n1==h) nxj=(long double)(*(*h).mp1).n1x;
			if ((*(*h).mp1).n2==h) nxj=(long double)(*(*h).mp1).n2x;
			xj=(*(*h).mp1).x-nxj*stl;
			yj=(*(*h).mp1).y-nxj*shl;
			lj=0.5L*(*(*h).st1).len;
			sg1=(*(*(*h).st1).rod).sit;
			cg1=(*(*(*h).st1).rod).cot;
			// Common for cases A,B,C and different for D
			if ((*h).st3)
			{
				// Cases A,B,C: at least 3 sticks
				sg2=(*(*(*h).st3).rod).sit;
				cg2=(*(*(*h).st3).rod).cot;
			}
			else
			{
				// Case D: 2 sticks on different rods
				sg2=(*(*(*h).st2).rod).sit;
				cg2=(*(*(*h).st2).rod).cot;
			}
			// Common for cases A and B, different for C and D
			if (((*h).st3)&&(((*(*h).st1).rod)==((*(*h).st2).rod)))
			{
				// Cases A,B: first two sticks on same rod
				if ((*(*h).mp2).n1==h) nxk=(long double)(*(*h).mp2).n1x;
				if ((*(*h).mp2).n2==h) nxk=(long double)(*(*h).mp2).n2x;
				if ((*(*h).mp3).n1==h) nxl=(long double)(*(*h).mp3).n1x;
				if ((*(*h).mp3).n2==h) nxl=(long double)(*(*h).mp3).n2x;
				xk=(*(*h).mp2).x-nxk*stl;
				yk=(*(*h).mp2).y-nxk*shl;
				xl=(*(*h).mp3).x-nxl*stl;
				yl=(*(*h).mp3).y-nxl*shl;
				lk=0.5L*(*(*h).st2).len;
				ll=0.5L*(*(*h).st3).len;
			}
			else
			{
				// Cases C,D: stick 1 and 2 on different rods
				if ((*(*h).mp2).n1==h) nxl=(long double)(*(*h).mp2).n1x;
				if ((*(*h).mp2).n2==h) nxl=(long double)(*(*h).mp2).n2x;
				xk=2.0L*(*h).x-xj;
				yk=2.0L*(*h).y-yj;
				xl=(*(*h).mp2).x-nxl*stl;
				yl=(*(*h).mp2).y-nxl*shl;
				lk=lj;
				ll=0.5L*(*(*h).st2).len;
			}
			// Only for case A, different for cases B,C,D
			if ((*h).st4)
			{
				// Case A: 4 sticks
				if ((*(*h).mp4).n1==h) nxm=(long double)(*(*h).mp4).n1x;
				if ((*(*h).mp4).n2==h) nxm=(long double)(*(*h).mp4).n2x;
				xm=(*(*h).mp4).x-nxm*stl;
				ym=(*(*h).mp4).y-nxm*shl;
				lm=0.5L*(*(*h).st4).len;
			}
			else
			{
				// Cases B,C,D: less than 4 sticks
				// Common for cases B and D, different for C
				if (((*h).st3) && (  ((*(*h).st2).rod)==((*(*h).st3).rod)  ))
				{
					// Case C: 3 sticks, sticks 2 and 3 on same rod
					if ((*(*h).mp3).n1==h) nxm=(long double)(*(*h).mp3).n1x;
					if ((*(*h).mp3).n2==h) nxm=(long double)(*(*h).mp3).n2x;
					xm=(*(*h).mp3).x-nxm*stl;
					ym=(*(*h).mp3).y-nxm*shl;
					lm=0.5L*(*(*h).st3).len;
				}
				else
				{
					// Cases B and D: either 3 sticks with sticks 1 and 2 on same rod,
					// or 2 sticks on different rods
					xm=2.0L*(*h).x-xl;
					ym=2.0L*(*h).y-yl;
					lm=ll;
				}
			}
			// Now we can compute the energy at this hinge
			// First we compute the angle change dg
			dg=0.0L;
			dg=dg+(*h).x*(sg2*(ll-lm)/(ll*lm)-sg1*(lj-lk)/(lj*lk));
			dg=dg+(*h).y*(cg2*(lm-ll)/(ll*lm)-cg1*(lk-lj)/(lj*lk));
			dg=dg-xj*sg1*lk/(lj*(lj+lk));
			dg=dg+yj*cg1*lk/(lj*(lj+lk));
			dg=dg+xk*sg1*lj/(lk*(lj+lk));
			dg=dg-yk*cg1*lj/(lk*(lj+lk));
			dg=dg+xl*sg2*lm/(ll*(ll+lm));
			dg=dg-yl*cg2*lm/(ll*(ll+lm));
			dg=dg-xm*sg2*ll/(lm*(ll+lm));
			dg=dg+ym*cg2*ll/(lm*(ll+lm));
			// Now the energy
			energy=energy+0.5L*GI*dg*dg;
		}
		h=(*h).next;
	}
	return energy;
}

// **************************

long double energyfroma(int hno, int mpno)
// Computes the energy from matrix Am, vector bv and the constant
{
	int qqx;
	struct hinges *hh;
	struct midpoints *mm;
	long double *vecvec;
	long double energiaabol;
	long int k; //NRsp 

	qqx=-1;
	vecvec=(long double *)calloc(2*(hno+mpno),sizeof(long double));
	hh=hingestart;
	while (hh)
	{
		qqx++;
		vecvec[qqx]=(*hh).x;
		qqx++;
		vecvec[qqx]=(*hh).y;
		hh=(*hh).next;
	}
	mm=midpointstart;
	while (mm)
	{
		qqx++;
		vecvec[qqx]=(*mm).x;
		qqx++;
		vecvec[qqx]=(*mm).y;
		mm=(*mm).next;
	}
	// Add the constant term first
	energiaabol=cconst;
	// Compute xAx/2
		for (k=1;k<=kspr;k++)
	{
		energiaabol=energiaabol+vecvec[Aspri[k]-1]*Aspr[k]*vecvec[Asprj[k]-1]/2.0L;
	}
	// Compute -bx
	for (qqx=0;qqx<2*(hno+mpno);qqx++)
	{
		energiaabol=energiaabol-bv[qqx]*vecvec[qqx];
	}
	return energiaabol;
}

// **************************

int conjgradvars(int hno, int mpno, long double shl, long double stl)
// Fills in the data of matrix A and vector b to solve Ax=b, then solution x
// gives minimizer of energy H(x)=xAx/2-bx+c.
{
	int S,SS,i,j,k,l,m,cont1,cont2,nix,nkx;
	long double eal0,eil0,sit,cot,ss,cc,sc,l1,l2,l1l2,DD;
	long double nxj,nxk,nxl,nxm,dxj,dyj,dxk,dyk,dxl,dyl,dxm,dym,lj,lk,ll,lm,sg1,cg1,sg2,cg2,xi,yi,xj,yj,xk,yk,xl,yl,xm,ym,om;
	struct midpoints *p;
	struct hinges *h;
	// Variables for spare matrix

	kspr=0;
	ksprMAX=0;
	ksprSTP=1024;
	// Total number of hinges and midpoints
	S=hno+mpno;
	SS=2*S;
	// Constant term
	cconst=0.0L;
	// Cycle along all midpoints to compile Am, bv and cconst from the attached half sticks
	// Contribution of stretch
	p=midpointstart;
	while (p)
	{
		eal0=0.5L*(*(*p).st).len;    // length of half-stick
		eal0=EA/eal0;               // EA/l0
		i=(*p).ser;                 // Serial number of actual midpoint
		sit=-(*(*(*p).st).rod).sit;  // Sine of rod angle
		cot=-(*(*(*p).st).rod).cot;  // Cosine of rod angle
		ss=eal0*sit*sit;
		cc=eal0*cot*cot;
		sc=eal0*sit*cot;
		// One end
		j=(*(*p).n1).ser;        // Serial number of neighboring hinge
		// Matrix Am
		mkAsp(SS,2*(SS*(hno+i-1)+hno+i-1),+cc);
		mkAsp(SS,2*(SS*(j-1)+j-1),+cc);
		mkAsp(SS,2*(SS*(hno+i)+hno+i)-1-SS,+ss);
		mkAsp(SS,2*(j*SS+j)-1-SS,+ss);
		mkAsp(SS,2*(SS*(hno+i-1)+j-1),-cc);
		mkAsp(SS,2*(SS*(j-1)+hno+i-1),-cc);
		mkAsp(SS,2*(SS*(hno+i)+j)-1-SS,-ss);
		mkAsp(SS,2*(j*SS+hno+i)-1-SS,-ss);
		mkAsp(SS,2*(SS*(hno+i-1)+hno+i)-1,+sc);
		mkAsp(SS,2*(SS*(hno+i)+hno+i-1)-SS,+sc);
		mkAsp(SS,2*(SS*(j-1)+j)-1,+sc);
		mkAsp(SS,2*(j*SS+j-1)-SS,+sc);
		mkAsp(SS,2*(SS*(j-1)+hno+i)-1,-sc);
		mkAsp(SS,2*(SS*(hno+i)+j-1)-SS,-sc);
		mkAsp(SS,2*(SS*(hno+i-1)+j)-1,-sc);
		mkAsp(SS,2*(j*SS+hno+i-1)-SS,-sc);
		// Vector bv
		DD=(cot*stl*((long double)((*p).n1x))+sit*shl*((long double)((*p).n1x)));
		bv[2*(hno+i-1)]=bv[2*(hno+i-1)]+eal0*cot*DD;
		bv[2*j-2]=bv[2*j-2]-eal0*cot*DD;
		bv[2*(hno+i)-1]=bv[2*(hno+i)-1]+eal0*sit*DD;
		bv[2*j-1]=bv[2*j-1]-eal0*sit*DD;
		// Constant c
		cconst=cconst+0.5L*eal0*DD*DD;
		// Other end
		// Here the sign of the angle can be changed if needed!
		// DO WE HAVE TO DO IT??? 
		// IT EFFECTIVELY CHANGES THE SIGN OF bv (?)
		// WE DO NOT DO IT NOW...
//		sit=-sit;
//		cot=-cot;
		j=(*(*p).n2).ser;        // Serial number of neighboring hinge
		// Matrix Am
		mkAsp(SS,2*(SS*(hno+i-1)+hno+i-1),+cc);
		mkAsp(SS,2*(SS*(j-1)+j-1),+cc);
		mkAsp(SS,2*(SS*(hno+i)+hno+i)-1-SS,+ss);
		mkAsp(SS,2*(j*SS+j)-1-SS,+ss);
		mkAsp(SS,2*(SS*(hno+i-1)+j-1),-cc);
		mkAsp(SS,2*(SS*(j-1)+hno+i-1),-cc);
		mkAsp(SS,2*(SS*(hno+i)+j)-1-SS,-ss);
		mkAsp(SS,2*(j*SS+hno+i)-1-SS,-ss);
		mkAsp(SS,2*(SS*(hno+i-1)+hno+i)-1,+sc);
		mkAsp(SS,2*(SS*(hno+i)+hno+i-1)-SS,+sc);
		mkAsp(SS,2*(SS*(j-1)+j)-1,+sc);
		mkAsp(SS,2*(j*SS+j-1)-SS,+sc);
		mkAsp(SS,2*(SS*(j-1)+hno+i)-1,-sc);
		mkAsp(SS,2*(SS*(hno+i)+j-1)-SS,-sc);
		mkAsp(SS,2*(SS*(hno+i-1)+j)-1,-sc);
		mkAsp(SS,2*(j*SS+hno+i-1)-SS,-sc);
		// Vector bv
		DD=(cot*stl*((long double)((*p).n2x))+sit*shl*((long double)((*p).n2x)));
		bv[2*(hno+i-1)]=bv[2*(hno+i-1)]+eal0*cot*DD;
		bv[2*j-2]=bv[2*j-2]-eal0*cot*DD;
		bv[2*(hno+i)-1]=bv[2*(hno+i)-1]+eal0*sit*DD;
		bv[2*j-1]=bv[2*j-1]-eal0*sit*DD;
		// Constant c
		cconst=cconst+0.5L*eal0*DD*DD;
		p=(*p).next;
	}
	// Contribution of bending
	// Cycle along all midpoints
	p=midpointstart;
	while (p)
	{
		j=(*p).ser;                 // Serial number of actual midpoint
		i=(*(*p).n1).ser;           // Serial number of one of the neighboring hinge
		k=(*(*p).n2).ser;           // Serial number of the other neighboring hinge
		sit=(*(*(*p).st).rod).sit;  // Sine of rod angle
		cot=(*(*(*p).st).rod).cot;  // Cosine of rod angle
		ss=sit*sit;
		cc=cot*cot;
		sc=sit*cot;
		l1=0.5L*(*(*p).st).len;      // Length of half stick
		l2=l1;
		l1l2=1.0L/l1+1.0L/l2;
		// Check if there is a stick on the other side 
		// of the neighbouring hinges: if not then energy is increased
		// so that later we need to consider only hinges with sticks
		// on both sides. We set eil0 (l) accordingly
		// Hinge on one end:
		cont1=1;
		if (((*(*p).n1).st1)==((*p).st))
		{
			if (((*(*p).n1).st2)==0)
			{
				cont1=0;
			}
			else
			{
				if (((*(*(*p).n1).st2).rod)!=((*(*p).st).rod))
				{
					cont1=0;
				}
			}
		}
		if (((*(*p).n1).st2)==((*p).st))
		{
			if (((*(*(*p).n1).st1).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n1).st3)==0)
				{
					cont1=0;
				}
				else
				{
					if (((*(*(*p).n1).st3).rod)!=((*(*p).st).rod))
					{
						printf("(A) ILYEN ESET NEM LEHET!!!\n");
						cont1=0;
					}
				}
			}
		}
		if (((*(*p).n1).st3)==((*p).st))
		{
			if (((*(*(*p).n1).st2).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n1).st4)==0)
				{
					cont1=0;
				}
				else
				{
					if (((*(*(*p).n1).st4).rod)!=((*(*p).st).rod))
					{
						printf("(B) ILYEN ESET NEM LEHET!!!\n");
						cont1=0;
					}
				}
			}
		}
		if (((*(*p).n1).st4)==((*p).st))
		{
			if (((*(*(*p).n1).st3).rod)!=((*(*p).st).rod))
			{
				cont1=0;
			}
		}
		// Hinge on other end:
		cont2=1;
		if (((*(*p).n2).st1)==((*p).st))
		{
			if (((*(*p).n2).st2)==0)
			{
				cont2=0;
			}
			else
			{
				if (((*(*(*p).n2).st2).rod)!=((*(*p).st).rod))
				{
					cont2=0;
				}
			}
		}
		if (((*(*p).n2).st2)==((*p).st))
		{
			if (((*(*(*p).n2).st1).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n2).st3)==0)
				{
					cont2=0;
				}
				else
				{
					if (((*(*(*p).n2).st3).rod)!=((*(*p).st).rod))
					{
						printf("(C) ILYEN ESET NEM LEHET!!!\n");
						cont2=0;
					}
				}
			}
		}
		if (((*(*p).n2).st3)==((*p).st))
		{
			if (((*(*(*p).n2).st2).rod)!=((*(*p).st).rod))
			{
				if (((*(*p).n2).st4)==0)
				{
					cont2=0;
				}
				else
				{
					if (((*(*(*p).n2).st4).rod)!=((*(*p).st).rod))
					{
						printf("(D) ILYEN ESET NEM LEHET!!!\n");
						cont2=0;
					}
				}
			}
		}
		if (((*(*p).n2).st4)==((*p).st))
		{
			if (((*(*(*p).n2).st3).rod)!=((*(*p).st).rod))
			{
				cont2=0;
			}
		}
		// Now we can compute eil0: l then kappa/l:
		eil0=l1*(0.5L+((long double)(cont1+cont2+cont1*cont2))/6.0L);
		eil0=EI/eil0;
		// Matrix Am
		mkAsp(SS,2*(SS*(i-1)+i-1),+eil0*ss/l1/l1);
		mkAsp(SS,2*((SS+1)*(hno+j-1)),+eil0*ss*l1l2*l1l2);
		mkAsp(SS,2*(SS*(k-1)+k-1),+eil0*ss/l2/l2);

		mkAsp(SS,2*(i*SS+i)-1-SS,+eil0*cc/l1/l1);
		mkAsp(SS,2*(SS*(hno+j)+hno+j)-1-SS,+eil0*cc*l1l2*l1l2);
		mkAsp(SS,2*(k*SS+k)-1-SS,+eil0*cc/l2/l2);

		mkAsp(SS,2*(SS*(i-1)+hno+j-1),-eil0*ss*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j-1)+i-1),-eil0*ss*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j-1)+k-1),-eil0*ss*l1l2/l2);
		mkAsp(SS,2*(SS*(k-1)+hno+j-1),-eil0*ss*l1l2/l2);

		mkAsp(SS,2*(SS*(i-1)+k-1),+eil0*ss/l1/l2);
		mkAsp(SS,2*(SS*(k-1)+i-1),+eil0*ss/l1/l2);

		mkAsp(SS,2*(i*SS+hno+j)-1-SS,-eil0*cc*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j)+i)-1-SS,-eil0*cc*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j)+k)-1-SS,-eil0*cc*l1l2/l2);
		mkAsp(SS,2*(k*SS+hno+j)-1-SS,-eil0*cc*l1l2/l2);

		mkAsp(SS,2*(i*SS+k)-1-SS,+eil0*cc/l1/l2);
		mkAsp(SS,2*(k*SS+i)-1-SS,+eil0*cc/l1/l2);

		mkAsp(SS,2*(SS*(i-1)+i)-1,-eil0*sit*cot/l1/l1);
		mkAsp(SS,2*(i*SS+i-1)-SS,-eil0*sit*cot/l1/l1);
		mkAsp(SS,2*(SS*(k-1)+k)-1,-eil0*sit*cot/l2/l2);
		mkAsp(SS,2*(k*SS+k-1)-SS,-eil0*sit*cot/l2/l2);

		mkAsp(SS,2*(SS*(k-1)+i)-1,-eil0*sit*cot/l1/l2);
		mkAsp(SS,2*(i*SS+k-1)-SS,-eil0*sit*cot/l1/l2);
		mkAsp(SS,2*(SS*(i-1)+k)-1,-eil0*sit*cot/l1/l2);
		mkAsp(SS,2*(k*SS+i-1)-SS,-eil0*sit*cot/l1/l2);

		mkAsp(SS,2*(SS*(hno+j-1)+hno+j)-1,-eil0*sit*cot*l1l2*l1l2);
		mkAsp(SS,2*(SS*(hno+j)+hno+j-1)-SS,-eil0*sit*cot*l1l2*l1l2);

		mkAsp(SS,2*(SS*(i-1)+hno+j)-1,+eil0*sit*cot*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j)+i-1)-SS,+eil0*sit*cot*l1l2/l1);
		mkAsp(SS,2*(SS*(hno+j-1)+i)-1,+eil0*sit*cot*l1l2/l1);
		mkAsp(SS,2*(i*SS+hno+j-1)-SS,+eil0*sit*cot*l1l2/l1);

		mkAsp(SS,2*(SS*(hno+j-1)+k)-1,+eil0*sit*cot*l1l2/l2);
		mkAsp(SS,2*(k*SS+hno+j-1)-SS,+eil0*sit*cot*l1l2/l2);
		mkAsp(SS,2*(SS*(k-1)+hno+j)-1,+eil0*sit*cot*l1l2/l2);
		mkAsp(SS,2*(SS*(hno+j)+k-1)-SS,+eil0*sit*cot*l1l2/l2);

		// Vector bv
		DD=((sit*stl-cot*shl)*(((long double)((*p).n1x))/l1+((long double)((*p).n2x))/l2));

		bv[2*i-2]=bv[2*i-2]-eil0*DD*sit/l1;
		bv[2*(hno+j-1)]=bv[2*(hno+j-1)]+eil0*DD*sit*l1l2;
		bv[2*k-2]=bv[2*k-2]-eil0*DD*sit/l2;

		bv[2*i-1]=bv[2*i-1]+eil0*DD*cot/l1;
		bv[2*(hno+j)-1]=bv[2*(hno+j)-1]-eil0*DD*cot*l1l2;
		bv[2*k-1]=bv[2*k-1]+eil0*DD*cot/l2;

		// Constant c
		cconst=cconst+0.5L*eil0*DD*DD;
		p=(*p).next;
	}
	// Cycle along all hinges
	h=hingestart;
	while (h)
	{
		// Compute bend of sticks at actual hinge
		// Need only to consider sticks if there is a stick on
		// both sides of the hinge, dangling ends are treated already

		// Check if there is any stick at all
		if (((*h).st1)&&((*h).st2))
		{
			// There are two sticks, check if they belong to same rod
			if (((*(*h).st1).rod)==((*(*h).st2).rod))
			{
				// OK, sticks 1 and 2 are on the same rod, compute Am matrix
				//
				j=(*h).ser;
				i=(*(*h).mp1).ser;
				k=(*(*h).mp2).ser;
				l1=0.5L*(*(*h).st1).len;
				l2=0.5L*(*(*h).st2).len;
				l1l2=1.0L/l1+1.0L/l2;
				eil0=2.0*EI/(l1+l2);
				sit=(*(*(*h).st1).rod).sit;
				cot=(*(*(*h).st1).rod).cot;
				ss=sit*sit;
				cc=cot*cot;
				sc=sit*cot;

				mkAsp(SS,2*(SS+1)*(hno+i-1),+eil0*ss/l1/l1);
				mkAsp(SS,2*(SS+1)*(j-1),+eil0*ss*l1l2*l1l2);
				mkAsp(SS,2*(SS+1)*(hno+k-1),+eil0*ss/l2/l2);

				mkAsp(SS,2*(SS*(hno+i)+hno+i)-1-SS,+eil0*cc/l1/l1);
				mkAsp(SS,2*(j*SS+j)-1-SS,+eil0*cc*l1l2*l1l2);
				mkAsp(SS,2*(SS*(hno+k)+hno+k)-1-SS,+eil0*cc/l2/l2);

				mkAsp(SS,2*(SS*(hno+i-1)+j-1),-eil0*ss*l1l2/l1);
				mkAsp(SS,2*(SS*(j-1)+hno+i-1),-eil0*ss*l1l2/l1);
				mkAsp(SS,2*(SS*(j-1)+hno+k-1),-eil0*ss*l1l2/l2);
				mkAsp(SS,2*(SS*(hno+k-1)+j-1),-eil0*ss*l1l2/l2);

				mkAsp(SS,2*(SS*(hno+i-1)+hno+k-1),+eil0*ss/l1/l2);
				mkAsp(SS,2*(SS*(hno+k-1)+hno+i-1),+eil0*ss/l1/l2);

				mkAsp(SS,2*(SS*(hno+i)+j)-1-SS,-eil0*cc*l1l2/l1);
				mkAsp(SS,2*(j*SS+hno+i)-1-SS,-eil0*cc*l1l2/l1);
				mkAsp(SS,2*(j*SS+hno+k)-1-SS,-eil0*cc*l1l2/l2);
				mkAsp(SS,2*(SS*(hno+k)+j)-1-SS,-eil0*cc*l1l2/l2);

				mkAsp(SS,2*(SS*(hno+i)+hno+k)-1-SS,+eil0*cc/l1/l2);
				mkAsp(SS,2*(SS*(hno+k)+hno+i)-1-SS,+eil0*cc/l1/l2);

				mkAsp(SS,2*(SS*(hno+i-1)+hno+i)-1,-eil0*sc/l1/l1);
				mkAsp(SS,2*(SS*(hno+i)+hno+i-1)-SS,-eil0*sc/l1/l1);
				mkAsp(SS,2*(SS*(hno+k-1)+hno+k)-1,-eil0*sc/l2/l2);
				mkAsp(SS,2*(SS*(hno+k)+hno+k-1)-SS,-eil0*sc/l2/l2);

				mkAsp(SS,2*(SS*(hno+k-1)+hno+i)-1,-eil0*sc/l1/l2);
				mkAsp(SS,2*(SS*(hno+i)+hno+k-1)-SS,-eil0*sc/l1/l2);
				mkAsp(SS,2*(SS*(hno+i-1)+hno+k)-1,-eil0*sc/l1/l2);
				mkAsp(SS,2*(SS*(hno+k)+hno+i-1)-SS,-eil0*sc/l1/l2);

				mkAsp(SS,2*(SS*(j-1)+j)-1,-eil0*sc*l1l2*l1l2);
				mkAsp(SS,2*(j*SS+j-1)-SS,-eil0*sc*l1l2*l1l2);

				mkAsp(SS,2*(SS*(hno+i-1)+j)-1,+eil0*sc*l1l2/l1);
				mkAsp(SS,2*(j*SS+hno+i-1)-SS,+eil0*sc*l1l2/l1);
				mkAsp(SS,2*(SS*(j-1)+hno+i)-1,+eil0*sc*l1l2/l1);
				mkAsp(SS,2*(SS*(hno+i)+j-1)-SS,+eil0*sc*l1l2/l1);

				mkAsp(SS,2*(SS*(j-1)+hno+k)-1,+eil0*sc*l1l2/l2);
				mkAsp(SS,2*(SS*(hno+k)+j-1)-SS,+eil0*sc*l1l2/l2);
				mkAsp(SS,2*(SS*(hno+k-1)+j)-1,+eil0*sc*l1l2/l2);
				mkAsp(SS,2*(j*SS+hno+k-1)-SS,+eil0*sc*l1l2/l2);

				// Vector bv
				if ((*(*h).mp1).n1==h)
				{
					nix=(*(*h).mp1).n1x;
				}
				else
				{
					nix=(*(*h).mp1).n2x;
				}
				if ((*(*h).mp2).n1==h)
				{
					nkx=(*(*h).mp2).n1x;
				}
				else
				{
					nkx=(*(*h).mp2).n2x;
				}
				DD=((-sit*stl+cot*shl)*(((long double)nix)/l1+((long double)nkx)/l2));

				bv[2*(hno+i-1)]=bv[2*(hno+i-1)]-eil0*DD*sit/l1;
				bv[2*j-2]=bv[2*j-2]+eil0*DD*sit*l1l2;
				bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-eil0*DD*sit/l2;

				bv[2*(hno+i)-1]=bv[2*(hno+i)-1]+eil0*DD*cot/l1;
				bv[2*j-1]=bv[2*j-1]-eil0*DD*cot*l1l2;
				bv[2*(hno+k)-1]=bv[2*(hno+k)-1]+eil0*DD*cot/l2;

				// Constant c
				cconst=cconst+0.5L*eil0*DD*DD;
		
				// Check if there are more sticks
				if (((*h).st3)&&((*h).st4))
				{
					// There are two more sticks, are they on the same rod?
					if (((*(*h).st3).rod)==((*(*h).st4).rod))
					{
						// OK, sticks 3 and 4 are on the same rod, compute Am matrix
						j=(*h).ser;
						i=(*(*h).mp3).ser;
						k=(*(*h).mp4).ser;
						l1=0.5L*(*(*h).st3).len;
						l2=0.5L*(*(*h).st4).len;
						l1l2=1.0L/l1+1.0L/l2;
						eil0=2.0L*EI/(l1+l2);
						sit=(*(*(*h).st3).rod).sit;
						cot=(*(*(*h).st3).rod).cot;
						ss=sit*sit;
						cc=cot*cot;
						sc=sit*cot;

						mkAsp(SS,2*(SS+1)*(hno+i-1),+eil0*ss/l1/l1);
						mkAsp(SS,2*(SS+1)*(j-1),+eil0*ss*l1l2*l1l2);
						mkAsp(SS,2*(SS+1)*(hno+k-1),+eil0*ss/l2/l2);

						mkAsp(SS,2*(SS*(hno+i)+hno+i)-1-SS,+eil0*cc/l1/l1);
						mkAsp(SS,2*(j*SS+j)-1-SS,+eil0*cc*l1l2*l1l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+k)-1-SS,+eil0*cc/l2/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+j-1),-eil0*ss*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+i-1),-eil0*ss*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+k-1),-eil0*ss*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+j-1),-eil0*ss*l1l2/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+hno+k-1),+eil0*ss/l1/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+hno+i-1),+eil0*ss/l1/l2);

						mkAsp(SS,2*(SS*(hno+i)+j)-1-SS,-eil0*cc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+i)-1-SS,-eil0*cc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+k)-1-SS,-eil0*cc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+j)-1-SS,-eil0*cc*l1l2/l2);

						mkAsp(SS,2*(SS*(hno+i)+hno+k)-1-SS,+eil0*cc/l1/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+i)-1-SS,+eil0*cc/l1/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+hno+i)-1,-eil0*sc/l1/l1);
						mkAsp(SS,2*(SS*(hno+i)+hno+i-1)-SS,-eil0*sc/l1/l1);
						mkAsp(SS,2*(SS*(hno+k-1)+hno+k)-1,-eil0*sc/l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+k-1)-SS,-eil0*sc/l2/l2);

						mkAsp(SS,2*(SS*(hno+k-1)+hno+i)-1,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+i)+hno+k-1)-SS,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+i-1)+hno+k)-1,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+i-1)-SS,-eil0*sc/l1/l2);

						mkAsp(SS,2*(SS*(j-1)+j)-1,-eil0*sc*l1l2*l1l2);
						mkAsp(SS,2*(j*SS+j-1)-SS,-eil0*sc*l1l2*l1l2);

						mkAsp(SS,2*(SS*(hno+i-1)+j)-1,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+i-1)-SS,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+i)-1,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(SS*(hno+i)+j-1)-SS,+eil0*sc*l1l2/l1);

						mkAsp(SS,2*(SS*(j-1)+hno+k)-1,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+j-1)-SS,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+j)-1,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(j*SS+hno+k-1)-SS,+eil0*sc*l1l2/l2);

						// Vector bv
						if ((*(*h).mp3).n1==h)
						{
							nix=(*(*h).mp3).n1x;
						}
						else
						{
							nix=(*(*h).mp3).n2x;
						}
						if ((*(*h).mp4).n1==h)
						{
							nkx=(*(*h).mp4).n1x;
						}
						else
						{
							nkx=(*(*h).mp4).n2x;
						}
						DD=((-sit*stl+cot*shl)*(((long double)nix)/l1+((long double)nkx)/l2));

						bv[2*(hno+i-1)]=bv[2*(hno+i-1)]-eil0*DD*sit/l1;
						bv[2*j-2]=bv[2*j-2]+eil0*DD*sit*l1l2;
						bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-eil0*DD*sit/l2;

						bv[2*(hno+i)-1]=bv[2*(hno+i)-1]+eil0*DD*cot/l1;
						bv[2*j-1]=bv[2*j-1]-eil0*DD*cot*l1l2;
						bv[2*(hno+k)-1]=bv[2*(hno+k)-1]+eil0*DD*cot/l2;

						// Constant c
						cconst=cconst+0.5L*eil0*DD*DD;
					}
				}
      		}
			else
			{
				// So first two sticks not on the same rod, are there more sticks?
				if ((*h).st3)
				{
					// Are sticks 2 and 3 on the same rod?
					if (((*(*h).st2).rod)==((*(*h).st3).rod))
					{
						// OK, sticks 2 and 3 are on the same rod, compute Am matrix
						j=(*h).ser;
						i=(*(*h).mp2).ser;
						k=(*(*h).mp3).ser;
						l1=0.5L*(*(*h).st2).len;
						l2=0.5L*(*(*h).st3).len;
						l1l2=1.0L/l1+1.0L/l2;
						eil0=2.0L*EI/(l1+l2);
						sit=(*(*(*h).st2).rod).sit;
						cot=(*(*(*h).st2).rod).cot;
						ss=sit*sit;
						cc=cot*cot;
						sc=sit*cot;

						mkAsp(SS,2*(SS+1)*(hno+i-1),+eil0*ss/l1/l1);
						mkAsp(SS,2*(SS+1)*(j-1),+eil0*ss*l1l2*l1l2);
						mkAsp(SS,2*(SS+1)*(hno+k-1),+eil0*ss/l2/l2);

						mkAsp(SS,2*(SS*(hno+i)+hno+i)-1-SS,+eil0*cc/l1/l1);
						mkAsp(SS,2*(j*SS+j)-1-SS,+eil0*cc*l1l2*l1l2);

						mkAsp(SS,2*(SS*(hno+k)+hno+k)-1-SS,+eil0*cc/l2/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+j-1),-eil0*ss*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+i-1),-eil0*ss*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+k-1),-eil0*ss*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+j-1),-eil0*ss*l1l2/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+hno+k-1),+eil0*ss/l1/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+hno+i-1),+eil0*ss/l1/l2);

						mkAsp(SS,2*(SS*(hno+i)+j)-1-SS,-eil0*cc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+i)-1-SS,-eil0*cc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+k)-1-SS,-eil0*cc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+j)-1-SS,-eil0*cc*l1l2/l2);

						mkAsp(SS,2*(SS*(hno+i)+hno+k)-1-SS,+eil0*cc/l1/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+i)-1-SS,+eil0*cc/l1/l2);

						mkAsp(SS,2*(SS*(hno+i-1)+hno+i)-1,-eil0*sc/l1/l1);
						mkAsp(SS,2*(SS*(hno+i)+hno+i-1)-SS,-eil0*sc/l1/l1);
						mkAsp(SS,2*(SS*(hno+k-1)+hno+k)-1,-eil0*sc/l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+k-1)-SS,-eil0*sc/l2/l2);

						mkAsp(SS,2*(SS*(hno+k-1)+hno+i)-1,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+i)+hno+k-1)-SS,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+i-1)+hno+k)-1,-eil0*sc/l1/l2);
						mkAsp(SS,2*(SS*(hno+k)+hno+i-1)-SS,-eil0*sc/l1/l2);

						mkAsp(SS,2*(SS*(j-1)+j)-1,-eil0*sc*l1l2*l1l2);
						mkAsp(SS,2*(j*SS+j-1)-SS,-eil0*sc*l1l2*l1l2);

						mkAsp(SS,2*(SS*(hno+i-1)+j)-1,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(j*SS+hno+i-1)-SS,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(SS*(j-1)+hno+i)-1,+eil0*sc*l1l2/l1);
						mkAsp(SS,2*(SS*(hno+i)+j-1)-SS,+eil0*sc*l1l2/l1);

						mkAsp(SS,2*(SS*(j-1)+hno+k)-1,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k)+j-1)-SS,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(SS*(hno+k-1)+j)-1,+eil0*sc*l1l2/l2);
						mkAsp(SS,2*(j*SS+hno+k-1)-SS,+eil0*sc*l1l2/l2);

						// Vector bv
						if ((*(*h).mp2).n1==h)
						{
							nix=(*(*h).mp2).n1x;
						}
						else
						{
							nix=(*(*h).mp2).n2x;
						}
						if ((*(*h).mp3).n1==h)
						{
							nkx=(*(*h).mp3).n1x;
						}
						else
						{
							nkx=(*(*h).mp3).n2x;
						}
						DD=((-sit*stl+cot*shl)*(((long double)nix)/l1+((long double)nkx)/l2));

						bv[2*(hno+i-1)]=bv[2*(hno+i-1)]-eil0*DD*sit/l1;
						bv[2*j-2]=bv[2*j-2]+eil0*DD*sit*l1l2;
						bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-eil0*DD*sit/l2;

						bv[2*(hno+i)-1]=bv[2*(hno+i)-1]+eil0*DD*cot/l1;
						bv[2*j-1]=bv[2*j-1]-eil0*DD*cot*l1l2;
						bv[2*(hno+k)-1]=bv[2*(hno+k)-1]+eil0*DD*cot/l2;

						// Constant c
						cconst=cconst+0.5L*eil0*DD*DD;
					}
				}

			}
		}
		h=(*h).next;
	}

	// Contribution of shear
	// Cycle along all hinges
	h=hingestart;
	while (h)
	{
		// Check how many sticks join here, need to have two connecting rods here with valid sticks
		if (((*h).st3) || (((*h).st2)&&(((*(*h).st1).rod)!=((*(*h).st2).rod))))
		{
			// We have at least 3 sticks or two sticks on different rods
			// Hence we have shear energy at this hinge
			// Common for all cases:
			i=(*h).ser;         // Serial number of actual hinge
			j=(*(*h).mp1).ser;  // Serial number of first midpoint
			if ((*(*h).mp1).n1==h) nxj=(long double)(*(*h).mp1).n1x;
			if ((*(*h).mp1).n2==h) nxj=(long double)(*(*h).mp1).n2x;
			dxj=-nxj*stl;
			dyj=-nxj*shl;
			lj=0.5L*(*(*h).st1).len;
			sg1=(*(*(*h).st1).rod).sit;
			cg1=(*(*(*h).st1).rod).cot;
			// Common for cases A,B,C and different for D
			if ((*h).st3)
			{
				// Cases A,B,C: at least 3 sticks
				sg2=(*(*(*h).st3).rod).sit;
				cg2=(*(*(*h).st3).rod).cot;
			}
			else
			{
				// Case D: 2 sticks on different rods
				sg2=(*(*(*h).st2).rod).sit;
				cg2=(*(*(*h).st2).rod).cot;
			}
			// Common for cases A and B, different for C and D
			if (((*h).st3)&&(((*(*h).st1).rod)==((*(*h).st2).rod)))
			{
				// Cases A,B: first two sticks on same rod
				k=(*(*h).mp2).ser;
				l=(*(*h).mp3).ser;
				if ((*(*h).mp2).n1==h) nxk=(long double)(*(*h).mp2).n1x;
				if ((*(*h).mp2).n2==h) nxk=(long double)(*(*h).mp2).n2x;
				if ((*(*h).mp3).n1==h) nxl=(long double)(*(*h).mp3).n1x;
				if ((*(*h).mp3).n2==h) nxl=(long double)(*(*h).mp3).n2x;
				dxk=-nxk*stl;
				dyk=-nxk*shl;
				dxl=-nxl*stl;
				dyl=-nxl*shl;
				lk=0.5L*(*(*h).st2).len;
				ll=0.5L*(*(*h).st3).len;
			}
			else
			{
				// Cases C,D: stick 1 and 2 on different rods
				l=(*(*h).mp2).ser;
				if ((*(*h).mp2).n1==h) nxl=(long double)(*(*h).mp2).n1x;
				if ((*(*h).mp2).n2==h) nxl=(long double)(*(*h).mp2).n2x;
				dxl=-nxl*stl;
				dyl=-nxl*shl;
				ll=0.5L*(*(*h).st2).len;
			}
			// Only for case A, different for cases B,C,D
			if ((*h).st4)
			{
				// Case A: 4 sticks
				m=(*(*h).mp4).ser;
				if ((*(*h).mp4).n1==h) nxm=(long double)(*(*h).mp4).n1x;
				if ((*(*h).mp4).n2==h) nxm=(long double)(*(*h).mp4).n2x;
				dxm=-nxm*stl;
				dym=-nxm*shl;
				lm=0.5L*(*(*h).st4).len;
			}
			else
			{
				// Cases B,C,D: less than 4 sticks
				// Common for cases B and D, different for C
				if (((*h).st3) && (((*(*h).st2).rod)==((*(*h).st3).rod)))
				{
					// Case C: 3 sticks, sticks 2 and 3 on same rod
					m=(*(*h).mp3).ser;
					if ((*(*h).mp3).n1==h) nxm=(long double)(*(*h).mp3).n1x;
					if ((*(*h).mp3).n2==h) nxm=(long double)(*(*h).mp3).n2x;
					dxm=-nxm*stl;
					dym=-nxm*shl;
					lm=0.5L*(*(*h).st3).len;
				}
				else
				{
					// Cases B and D: either 3 sticks with sticks 1 and 2 on same rod,
					// or 2 sticks on different rods
					// Nothing to do.
				}
			}
			// Compute now matrix A, vector b and constants c for each case
			// Case A
			if ((*h).st4)
			{
				// Case A: 4 sticks
				xi= sg2*(ll-lm)/(ll*lm)-sg1*(lj-lk)/(lj*lk);
				yi= cg2*(lm-ll)/(ll*lm)-cg1*(lk-lj)/(lj*lk);
				xj=-sg1*lk/(lj*(lj+lk));
				yj= cg1*lk/(lj*(lj+lk));
				xk= sg1*lj/(lk*(lj+lk));
				yk=-cg1*lj/(lk*(lj+lk));
				xl= sg2*lm/(ll*(ll+lm));
				yl=-cg2*lm/(ll*(ll+lm));
				xm=-sg2*ll/(lm*(ll+lm));
				ym= cg2*ll/(lm*(ll+lm));
				om=(xj*dxj+yj*dyj+xk*dxk+yk*dyk+xl*dxl+yl*dyl+xm*dxm+ym*dym);
				// Hinge-hinge types
				mkAsp(SS,2*(SS*(i-1)+i-1),+GI*xi*xi); // Xi Xi
				mkAsp(SS,2*(i*SS+i)-1-SS,+GI*yi*yi);  // Yi Yi
				mkAsp(SS,2*(SS*(i-1)+i)-1,+GI*xi*yi);  // Xi Yi
				mkAsp(SS,2*(i*SS+i-1)-SS,+GI*yi*xi);  // Yi Xi
				// Hinge-midpoint types
				mkAsp(SS,2*(SS*(i-1)+hno+j-1),+GI*xi*xj);  // Xi Xj
				mkAsp(SS,2*(SS*(i-1)+hno+k-1),+GI*xi*xk);  // Xi Xk
				mkAsp(SS,2*(SS*(i-1)+hno+l-1),+GI*xi*xl);  // Xi Xl
				mkAsp(SS,2*(SS*(i-1)+hno+m-1),+GI*xi*xm);  // Xi Xm

				mkAsp(SS,2*(i*SS+hno+j)-1-SS,+GI*yi*yj);  // Yi Yj
				mkAsp(SS,2*(i*SS+hno+k)-1-SS,+GI*yi*yk);  // Yi Yk
				mkAsp(SS,2*(i*SS+hno+l)-1-SS,+GI*yi*yl);  // Yi Yl
				mkAsp(SS,2*(i*SS+hno+m)-1-SS,+GI*yi*ym);  // Yi Ym

				mkAsp(SS,2*(SS*(i-1)+hno+j)-1,+GI*xi*yj);  // Xi Yj
				mkAsp(SS,2*(SS*(i-1)+hno+k)-1,+GI*xi*yk);  // Xi Yk
				mkAsp(SS,2*(SS*(i-1)+hno+l)-1,+GI*xi*yl);  // Xi Yl
				mkAsp(SS,2*(SS*(i-1)+hno+m)-1,+GI*xi*ym);  // Xi Ym

				mkAsp(SS,2*(i*SS+hno+j-1)-SS,+GI*yi*xj);  // Yi Xj
				mkAsp(SS,2*(i*SS+hno+k-1)-SS,+GI*yi*xk);  // Yi Xk
				mkAsp(SS,2*(i*SS+hno+l-1)-SS,+GI*yi*xl);  // Yi Xl
				mkAsp(SS,2*(i*SS+hno+m-1)-SS,+GI*yi*xm);  // Yi Xm

				// Midpoint-hinge types
				mkAsp(SS,2*(SS*(hno+j-1)+i-1),+GI*xj*xi);  // Xj Xi
				mkAsp(SS,2*(SS*(hno+k-1)+i-1),+GI*xk*xi);  // Xk Xi
				mkAsp(SS,2*(SS*(hno+l-1)+i-1),+GI*xl*xi);  // Xl Xi
				mkAsp(SS,2*(SS*(hno+m-1)+i-1),+GI*xm*xi);  // Xm Xi

				mkAsp(SS,2*(SS*(hno+j)+i)-1-SS,+GI*yj*yi);  // Yj Yi
				mkAsp(SS,2*(SS*(hno+k)+i)-1-SS,+GI*yk*yi);  // Yk Yi
				mkAsp(SS,2*(SS*(hno+l)+i)-1-SS,+GI*yl*yi);  // Yl Yi
				mkAsp(SS,2*(SS*(hno+m)+i)-1-SS,+GI*ym*yi);  // Ym Yi
				mkAsp(SS,2*(SS*(hno+j)+i-1)-SS,+GI*yj*xi);  // Yj Xi
				mkAsp(SS,2*(SS*(hno+k)+i-1)-SS,+GI*yk*xi);  // Yk Xi
				mkAsp(SS,2*(SS*(hno+l)+i-1)-SS,+GI*yl*xi);  // Yl Xi
				mkAsp(SS,2*(SS*(hno+m)+i-1)-SS,+GI*ym*xi);  // Ym Xi

				mkAsp(SS,2*(SS*(hno+j-1)+i)-1,+GI*xj*yi);  // Xj Yi
				mkAsp(SS,2*(SS*(hno+k-1)+i)-1,+GI*xk*yi);  // Xk Yi
				mkAsp(SS,2*(SS*(hno+l-1)+i)-1,+GI*xl*yi);  // Xl Yi
				mkAsp(SS,2*(SS*(hno+m-1)+i)-1,+GI*xm*yi);  // Xm Yi

				// Midpoint-midpoint types
				mkAsp(SS,2*(SS*(hno+j-1)+hno+j-1),+GI*xj*xj);  // Xj Xj
				mkAsp(SS,2*(SS*(hno+j-1)+hno+k-1),+GI*xj*xk);  // Xj Xk
				mkAsp(SS,2*(SS*(hno+j-1)+hno+l-1),+GI*xj*xl);  // Xj Xl
				mkAsp(SS,2*(SS*(hno+j-1)+hno+m-1),+GI*xj*xm);  // Xj Xm
				mkAsp(SS,2*(SS*(hno+j-1)+hno+j)-1,+GI*xj*yj);  // Xj Yj
				mkAsp(SS,2*(SS*(hno+j-1)+hno+k)-1,+GI*xj*yk);  // Xj Yk
				mkAsp(SS,2*(SS*(hno+j-1)+hno+l)-1,+GI*xj*yl);  // Xj Yl
				mkAsp(SS,2*(SS*(hno+j-1)+hno+m)-1,+GI*xj*ym);  // Xj Ym

				mkAsp(SS,2*(SS*(hno+k-1)+hno+j-1),+GI*xk*xj);  // Xk Xj
				mkAsp(SS,2*(SS*(hno+k-1)+hno+k-1),+GI*xk*xk);  // Xk Xk
				mkAsp(SS,2*(SS*(hno+k-1)+hno+l-1),+GI*xk*xl);  // Xk Xl
				mkAsp(SS,2*(SS*(hno+k-1)+hno+m-1),+GI*xk*xm);  // Xk Xm
				mkAsp(SS,2*(SS*(hno+k-1)+hno+j)-1,+GI*xk*yj);  // Xk Yj
				mkAsp(SS,2*(SS*(hno+k-1)+hno+k)-1,+GI*xk*yk);  // Xk Yk
				mkAsp(SS,2*(SS*(hno+k-1)+hno+l)-1,+GI*xk*yl);  // Xk Yl
				mkAsp(SS,2*(SS*(hno+k-1)+hno+m)-1,+GI*xk*ym);  // Xk Ym

				mkAsp(SS,2*(SS*(hno+l-1)+hno+j-1),+GI*xl*xj);  // Xl Xj
				mkAsp(SS,2*(SS*(hno+l-1)+hno+k-1),+GI*xl*xk);  // Xl Xk
				mkAsp(SS,2*(SS*(hno+l-1)+hno+l-1),+GI*xl*xl);  // Xl Xl
				mkAsp(SS,2*(SS*(hno+l-1)+hno+m-1),+GI*xl*xm);  // Xl Xm
				mkAsp(SS,2*(SS*(hno+l-1)+hno+j)-1,+GI*xl*yj);  // Xl Yj
				mkAsp(SS,2*(SS*(hno+l-1)+hno+k)-1,+GI*xl*yk);  // Xl Yk
				mkAsp(SS,2*(SS*(hno+l-1)+hno+l)-1,+GI*xl*yl);  // Xl Yl
				mkAsp(SS,2*(SS*(hno+l-1)+hno+m)-1,+GI*xl*ym);  // Xl Ym

				mkAsp(SS,2*(SS*(hno+m-1)+hno+j-1),+GI*xm*xj);  // Xm Xj
				mkAsp(SS,2*(SS*(hno+m-1)+hno+k-1),+GI*xm*xk);  // Xm Xk
				mkAsp(SS,2*(SS*(hno+m-1)+hno+l-1),+GI*xm*xl);  // Xm Xl
				mkAsp(SS,2*(SS*(hno+m-1)+hno+m-1),+GI*xm*xm);  // Xm Xm
				mkAsp(SS,2*(SS*(hno+m-1)+hno+j)-1,+GI*xm*yj);  // Xm Yj
				mkAsp(SS,2*(SS*(hno+m-1)+hno+k)-1,+GI*xm*yk);  // Xm Yk
				mkAsp(SS,2*(SS*(hno+m-1)+hno+l)-1,+GI*xm*yl);  // Xm Yl
				mkAsp(SS,2*(SS*(hno+m-1)+hno+m)-1,+GI*xm*ym);  // Xm Ym

				mkAsp(SS,2*(SS*(hno+j)+hno+j-1)-SS,+GI*yj*xj);  // Yj Xj
				mkAsp(SS,2*(SS*(hno+j)+hno+k-1)-SS,+GI*yj*xk);  // Yj Xk
				mkAsp(SS,2*(SS*(hno+j)+hno+l-1)-SS,+GI*yj*xl);  // Yj Xl
				mkAsp(SS,2*(SS*(hno+j)+hno+m-1)-SS,+GI*yj*xm);  // Yj Xm
				mkAsp(SS,2*(SS*(hno+j)+hno+j)-1-SS,+GI*yj*yj);  // Yj Yj
				mkAsp(SS,2*(SS*(hno+j)+hno+k)-1-SS,+GI*yj*yk);  // Yj Yk
				mkAsp(SS,2*(SS*(hno+j)+hno+l)-1-SS,+GI*yj*yl);  // Yj Yl
				mkAsp(SS,2*(SS*(hno+j)+hno+m)-1-SS,+GI*yj*ym);  // Yj Ym
                                                                
				mkAsp(SS,2*(SS*(hno+k)+hno+j-1)-SS,+GI*yk*xj);  // Yk Xj
				mkAsp(SS,2*(SS*(hno+k)+hno+k-1)-SS,+GI*yk*xk);  // Yk Xk
				mkAsp(SS,2*(SS*(hno+k)+hno+l-1)-SS,+GI*yk*xl);  // Yk Xl
				mkAsp(SS,2*(SS*(hno+k)+hno+m-1)-SS,+GI*yk*xm);  // Yk Xm
				mkAsp(SS,2*(SS*(hno+k)+hno+j)-1-SS,+GI*yk*yj);  // Yk Yj
				mkAsp(SS,2*(SS*(hno+k)+hno+k)-1-SS,+GI*yk*yk);  // Yk Yk
				mkAsp(SS,2*(SS*(hno+k)+hno+l)-1-SS,+GI*yk*yl);  // Yk Yl
				mkAsp(SS,2*(SS*(hno+k)+hno+m)-1-SS,+GI*yk*ym);  // Yk Ym
                                                                
				mkAsp(SS,2*(SS*(hno+l)+hno+j-1)-SS,+GI*yl*xj);  // Yl Xj
				mkAsp(SS,2*(SS*(hno+l)+hno+k-1)-SS,+GI*yl*xk);  // Yl Xk
				mkAsp(SS,2*(SS*(hno+l)+hno+l-1)-SS,+GI*yl*xl);  // Yl Xl
				mkAsp(SS,2*(SS*(hno+l)+hno+m-1)-SS,+GI*yl*xm);  // Yl Xm
				mkAsp(SS,2*(SS*(hno+l)+hno+j)-1-SS,+GI*yl*yj);  // Yl Yj
				mkAsp(SS,2*(SS*(hno+l)+hno+k)-1-SS,+GI*yl*yk);  // Yl Yk
				mkAsp(SS,2*(SS*(hno+l)+hno+l)-1-SS,+GI*yl*yl);  // Yl Yl
				mkAsp(SS,2*(SS*(hno+l)+hno+m)-1-SS,+GI*yl*ym);  // Yl Ym
                                                                
				mkAsp(SS,2*(SS*(hno+m)+hno+j-1)-SS,+GI*ym*xj);  // Ym Xj
				mkAsp(SS,2*(SS*(hno+m)+hno+k-1)-SS,+GI*ym*xk);  // Ym Xk
				mkAsp(SS,2*(SS*(hno+m)+hno+l-1)-SS,+GI*ym*xl);  // Ym Xl
				mkAsp(SS,2*(SS*(hno+m)+hno+m-1)-SS,+GI*ym*xm);  // Ym Xm
				mkAsp(SS,2*(SS*(hno+m)+hno+j)-1-SS,+GI*ym*yj);  // Ym Yj
				mkAsp(SS,2*(SS*(hno+m)+hno+k)-1-SS,+GI*ym*yk);  // Ym Yk
				mkAsp(SS,2*(SS*(hno+m)+hno+l)-1-SS,+GI*ym*yl);  // Ym Yl
				mkAsp(SS,2*(SS*(hno+m)+hno+m)-1-SS,+GI*ym*ym);  // Ym Ym

				bv[2*i-2]=bv[2*i-2]-GI*om*xi;  // Xi
				bv[2*i-1]=bv[2*i-1]-GI*om*yi;  // Yi

				bv[2*(hno+j-1)]=bv[2*(hno+j-1)]-GI*om*xj;  // Xj
				bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-GI*om*xk;  // Xk
				bv[2*(hno+l-1)]=bv[2*(hno+l-1)]-GI*om*xl;  // Xl
				bv[2*(hno+m-1)]=bv[2*(hno+m-1)]-GI*om*xm;  // Xm

				bv[2*(hno+j)-1]=bv[2*(hno+j)-1]-GI*om*yj;  // Yj
				bv[2*(hno+k)-1]=bv[2*(hno+k)-1]-GI*om*yk;  // Yk
				bv[2*(hno+l)-1]=bv[2*(hno+l)-1]-GI*om*yl;  // Yl
				bv[2*(hno+m)-1]=bv[2*(hno+m)-1]-GI*om*ym;  // Ym

				cconst=cconst+0.5L*GI*om*om;
			}
			else
			{
				// Case B
				if (((*h).st3)&&(((*(*h).st1).rod)==((*(*h).st2).rod)))
				{
					// Case B: 3 sticks, st1 and st2 on same rod
					xi=-sg2/ll-sg1*(lj-lk)/(lj*lk);
					yi= cg2/ll-cg1*(lk-lj)/(lj*lk);
					xj=-sg1*lk/(lj*(lj+lk));
					yj= cg1*lk/(lj*(lj+lk));
					xk= sg1*lj/(lk*(lj+lk));
					yk=-cg1*lj/(lk*(lj+lk));
					xl= sg2/ll;
					yl=-cg2/ll;
					om=(xj*dxj+yj*dyj+xk*dxk+yk*dyk+xl*dxl+yl*dyl);

					// Hinge-hinge types
					mkAsp(SS,2*(SS*(i-1)+i-1),+GI*xi*xi); // Xi Xi
					mkAsp(SS,2*(i*SS+i)-1-SS,+GI*yi*yi);  // Yi Yi
					mkAsp(SS,2*(SS*(i-1)+i)-1,+GI*xi*yi);  // Xi Yi
					mkAsp(SS,2*(i*SS+i-1)-SS,+GI*yi*xi);  // Yi Xi
	
					// Hinge-midpoint types
					mkAsp(SS,2*(SS*(i-1)+hno+j-1),+GI*xi*xj);  // Xi Xj
					mkAsp(SS,2*(SS*(i-1)+hno+k-1),+GI*xi*xk);  // Xi Xk
					mkAsp(SS,2*(SS*(i-1)+hno+l-1),+GI*xi*xl);  // Xi Xl
	
					mkAsp(SS,2*(i*SS+hno+j)-1-SS,+GI*yi*yj);  // Yi Yj
					mkAsp(SS,2*(i*SS+hno+k)-1-SS,+GI*yi*yk);  // Yi Yk
					mkAsp(SS,2*(i*SS+hno+l)-1-SS,+GI*yi*yl);  // Yi Yl
	
					mkAsp(SS,2*(SS*(i-1)+hno+j)-1,+GI*xi*yj);  // Xi Yj
					mkAsp(SS,2*(SS*(i-1)+hno+k)-1,+GI*xi*yk);  // Xi Yk
					mkAsp(SS,2*(SS*(i-1)+hno+l)-1,+GI*xi*yl);  // Xi Yl
	
					mkAsp(SS,2*(i*SS+hno+j-1)-SS,+GI*yi*xj);  // Yi Xj
					mkAsp(SS,2*(i*SS+hno+k-1)-SS,+GI*yi*xk);  // Yi Xk
					mkAsp(SS,2*(i*SS+hno+l-1)-SS,+GI*yi*xl);  // Yi Xl
	
					// Midpoint-hinge types
					mkAsp(SS,2*(SS*(hno+j-1)+i-1),+GI*xj*xi);  // Xj Xi
					mkAsp(SS,2*(SS*(hno+k-1)+i-1),+GI*xk*xi);  // Xk Xi
					mkAsp(SS,2*(SS*(hno+l-1)+i-1),+GI*xl*xi);  // Xl Xi
	
					mkAsp(SS,2*(SS*(hno+j)+i)-1-SS,+GI*yj*yi);  // Yj Yi
					mkAsp(SS,2*(SS*(hno+k)+i)-1-SS,+GI*yk*yi);  // Yk Yi
					mkAsp(SS,2*(SS*(hno+l)+i)-1-SS,+GI*yl*yi);  // Yl Yi
	
					mkAsp(SS,2*(SS*(hno+j)+i-1)-SS,+GI*yj*xi);  // Yj Xi
					mkAsp(SS,2*(SS*(hno+k)+i-1)-SS,+GI*yk*xi);  // Yk Xi
					mkAsp(SS,2*(SS*(hno+l)+i-1)-SS,+GI*yl*xi);  // Yl Xi
	
					mkAsp(SS,2*(SS*(hno+j-1)+i)-1,+GI*xj*yi);  // Xj Yi
					mkAsp(SS,2*(SS*(hno+k-1)+i)-1,+GI*xk*yi);  // Xk Yi
					mkAsp(SS,2*(SS*(hno+l-1)+i)-1,+GI*xl*yi);  // Xl Yi
	
					// Midpoint-midpoint types
					mkAsp(SS,2*(SS*(hno+j-1)+hno+j-1),+GI*xj*xj);  // Xj Xj
					mkAsp(SS,2*(SS*(hno+j-1)+hno+k-1),+GI*xj*xk);  // Xj Xk
					mkAsp(SS,2*(SS*(hno+j-1)+hno+l-1),+GI*xj*xl);  // Xj Xl
					mkAsp(SS,2*(SS*(hno+j-1)+hno+j)-1,+GI*xj*yj);  // Xj Yj
					mkAsp(SS,2*(SS*(hno+j-1)+hno+k)-1,+GI*xj*yk);  // Xj Yk
					mkAsp(SS,2*(SS*(hno+j-1)+hno+l)-1,+GI*xj*yl);  // Xj Yl
	
					mkAsp(SS,2*(SS*(hno+k-1)+hno+j-1),+GI*xk*xj);  // Xk Xj
					mkAsp(SS,2*(SS*(hno+k-1)+hno+k-1),+GI*xk*xk);  // Xk Xk
					mkAsp(SS,2*(SS*(hno+k-1)+hno+l-1),+GI*xk*xl);  // Xk Xl
					mkAsp(SS,2*(SS*(hno+k-1)+hno+j)-1,+GI*xk*yj);  // Xk Yj
					mkAsp(SS,2*(SS*(hno+k-1)+hno+k)-1,+GI*xk*yk);  // Xk Yk
					mkAsp(SS,2*(SS*(hno+k-1)+hno+l)-1,+GI*xk*yl);  // Xk Yl
	
					mkAsp(SS,2*(SS*(hno+l-1)+hno+j-1),+GI*xl*xj);  // Xl Xj
					mkAsp(SS,2*(SS*(hno+l-1)+hno+k-1),+GI*xl*xk);  // Xl Xk
					mkAsp(SS,2*(SS*(hno+l-1)+hno+l-1),+GI*xl*xl);  // Xl Xl
					mkAsp(SS,2*(SS*(hno+l-1)+hno+j)-1,+GI*xl*yj);  // Xl Yj
					mkAsp(SS,2*(SS*(hno+l-1)+hno+k)-1,+GI*xl*yk);  // Xl Yk
					mkAsp(SS,2*(SS*(hno+l-1)+hno+l)-1,+GI*xl*yl);  // Xl Yl

					mkAsp(SS,2*(SS*(hno+j)+hno+j-1)-SS,+GI*yj*xj);  // Yj Xj
					mkAsp(SS,2*(SS*(hno+j)+hno+k-1)-SS,+GI*yj*xk);  // Yj Xk
					mkAsp(SS,2*(SS*(hno+j)+hno+l-1)-SS,+GI*yj*xl);  // Yj Xl
					mkAsp(SS,2*(SS*(hno+j)+hno+j)-1-SS,+GI*yj*yj);  // Yj Yj
					mkAsp(SS,2*(SS*(hno+j)+hno+k)-1-SS,+GI*yj*yk);  // Yj Yk
					mkAsp(SS,2*(SS*(hno+j)+hno+l)-1-SS,+GI*yj*yl);  // Yj Yl
	
					mkAsp(SS,2*(SS*(hno+k)+hno+j-1)-SS,+GI*yk*xj);  // Yk Xj
					mkAsp(SS,2*(SS*(hno+k)+hno+k-1)-SS,+GI*yk*xk);  // Yk Xk
					mkAsp(SS,2*(SS*(hno+k)+hno+l-1)-SS,+GI*yk*xl);  // Yk Xl
					mkAsp(SS,2*(SS*(hno+k)+hno+j)-1-SS,+GI*yk*yj);  // Yk Yj
					mkAsp(SS,2*(SS*(hno+k)+hno+k)-1-SS,+GI*yk*yk);  // Yk Yk
					mkAsp(SS,2*(SS*(hno+k)+hno+l)-1-SS,+GI*yk*yl);  // Yk Yl
	
					mkAsp(SS,2*(SS*(hno+l)+hno+j-1)-SS,+GI*yl*xj);  // Yl Xj
					mkAsp(SS,2*(SS*(hno+l)+hno+k-1)-SS,+GI*yl*xk);  // Yl Xk
					mkAsp(SS,2*(SS*(hno+l)+hno+l-1)-SS,+GI*yl*xl);  // Yl Xl
					mkAsp(SS,2*(SS*(hno+l)+hno+j)-1-SS,+GI*yl*yj);  // Yl Yj
					mkAsp(SS,2*(SS*(hno+l)+hno+k)-1-SS,+GI*yl*yk);  // Yl Yk
					mkAsp(SS,2*(SS*(hno+l)+hno+l)-1-SS,+GI*yl*yl);  // Yl Yl
	
					bv[2*i-2]=bv[2*i-2]-GI*om*xi;  // Xi
					bv[2*i-1]=bv[2*i-1]-GI*om*yi;  // Yi
	
					bv[2*(hno+j-1)]=bv[2*(hno+j-1)]-GI*om*xj;  // Xj
					bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-GI*om*xk;  // Xk
					bv[2*(hno+l-1)]=bv[2*(hno+l-1)]-GI*om*xl;  // Xl
					// bv[2*(hno+m-1)]=bv[2*(hno+m-1)]-GI*om*xm;  // Xm
	
					bv[2*(hno+j)-1]=bv[2*(hno+j)-1]-GI*om*yj;  // Yj
					bv[2*(hno+k)-1]=bv[2*(hno+k)-1]-GI*om*yk;  // Yk
					bv[2*(hno+l)-1]=bv[2*(hno+l)-1]-GI*om*yl;  // Yl
					// bv[2*(hno+m)-1]=bv[2*(hno+m)-1]-GI*om*ym;  // Ym
	
					cconst=cconst+0.5L*GI*om*om;
				}
				else
				{
					// Case C
					if (((*h).st3)&&(((*(*h).st2).rod)==((*(*h).st3).rod)))
					{
						// Case C: 3 sticks, st2 and st3 on same rod
						xi= sg2*(ll-lm)/(ll*lm)+sg1/lj;
						yi= cg2*(lm-ll)/(ll*lm)-cg1/lj;
						xj=-sg1/lj;
						yj= cg1/lj;
						xl= sg2*lm/(ll*(ll+lm));
						yl=-cg2*lm/(ll*(ll+lm));
						xm=-sg2*ll/(lm*(ll+lm));
						ym= cg2*ll/(lm*(ll+lm));
						om=(xj*dxj+yj*dyj+xl*dxl+yl*dyl+xm*dxm+ym*dym);
						// Hinge-hinge types
						mkAsp(SS,2*(SS*(i-1)+i-1),+GI*xi*xi); // Xi Xi
						mkAsp(SS,2*(i*SS+i)-1-SS,+GI*yi*yi);  // Yi Yi
						mkAsp(SS,2*(SS*(i-1)+i)-1,+GI*xi*yi);  // Xi Yi
						mkAsp(SS,2*(i*SS+i-1)-SS,+GI*yi*xi);  // Yi Xi
		
						// Hinge-midpoint types
						mkAsp(SS,2*(SS*(i-1)+hno+j-1),+GI*xi*xj);  // Xi Xj
						mkAsp(SS,2*(SS*(i-1)+hno+l-1),+GI*xi*xl);  // Xi Xl
						mkAsp(SS,2*(SS*(i-1)+hno+m-1),+GI*xi*xm);  // Xi Xm
		
						mkAsp(SS,2*(i*SS+hno+j)-1-SS,+GI*yi*yj);  // Yi Yj
						mkAsp(SS,2*(i*SS+hno+l)-1-SS,+GI*yi*yl);  // Yi Yl
						mkAsp(SS,2*(i*SS+hno+m)-1-SS,+GI*yi*ym);  // Yi Ym
		
						mkAsp(SS,2*(SS*(i-1)+hno+j)-1,+GI*xi*yj);  // Xi Yj
						mkAsp(SS,2*(SS*(i-1)+hno+l)-1,+GI*xi*yl);  // Xi Yl
						mkAsp(SS,2*(SS*(i-1)+hno+m)-1,+GI*xi*ym);  // Xi Ym
		
						mkAsp(SS,2*(i*SS+hno+j-1)-SS,+GI*yi*xj);  // Yi Xj
						mkAsp(SS,2*(i*SS+hno+l-1)-SS,+GI*yi*xl);  // Yi Xl
						mkAsp(SS,2*(i*SS+hno+m-1)-SS,+GI*yi*xm);  // Yi Xm
		
						// Midpoint-hinge types
						mkAsp(SS,2*(SS*(hno+j-1)+i-1),+GI*xj*xi);  // Xj Xi
						mkAsp(SS,2*(SS*(hno+l-1)+i-1),+GI*xl*xi);  // Xl Xi
						mkAsp(SS,2*(SS*(hno+m-1)+i-1),+GI*xm*xi);  // Xm Xi
		
						mkAsp(SS,2*(SS*(hno+j)+i)-1-SS,+GI*yj*yi);  // Yj Yi
						mkAsp(SS,2*(SS*(hno+l)+i)-1-SS,+GI*yl*yi);  // Yl Yi
						mkAsp(SS,2*(SS*(hno+m)+i)-1-SS,+GI*ym*yi);  // Ym Yi
		
						mkAsp(SS,2*(SS*(hno+j)+i-1)-SS,+GI*yj*xi);  // Yj Xi
						mkAsp(SS,2*(SS*(hno+l)+i-1)-SS,+GI*yl*xi);  // Yl Xi
						mkAsp(SS,2*(SS*(hno+m)+i-1)-SS,+GI*ym*xi);  // Ym Xi
		
						mkAsp(SS,2*(SS*(hno+j-1)+i)-1,+GI*xj*yi);  // Xj Yi
						mkAsp(SS,2*(SS*(hno+l-1)+i)-1,+GI*xl*yi);  // Xl Yi
						mkAsp(SS,2*(SS*(hno+m-1)+i)-1,+GI*xm*yi);  // Xm Yi
		
						// Midpoint-midpoint types
						mkAsp(SS,2*(SS*(hno+j-1)+hno+j-1),+GI*xj*xj);  // Xj Xj
						mkAsp(SS,2*(SS*(hno+j-1)+hno+l-1),+GI*xj*xl);  // Xj Xl
						mkAsp(SS,2*(SS*(hno+j-1)+hno+m-1),+GI*xj*xm);  // Xj Xm
						mkAsp(SS,2*(SS*(hno+j-1)+hno+j)-1,+GI*xj*yj);  // Xj Yj
						mkAsp(SS,2*(SS*(hno+j-1)+hno+l)-1,+GI*xj*yl);  // Xj Yl
						mkAsp(SS,2*(SS*(hno+j-1)+hno+m)-1,+GI*xj*ym);  // Xj Ym
		
						mkAsp(SS,2*(SS*(hno+l-1)+hno+j-1),+GI*xl*xj);  // Xl Xj
						mkAsp(SS,2*(SS*(hno+l-1)+hno+l-1),+GI*xl*xl);  // Xl Xl
						mkAsp(SS,2*(SS*(hno+l-1)+hno+m-1),+GI*xl*xm);  // Xl Xm
						mkAsp(SS,2*(SS*(hno+l-1)+hno+j)-1,+GI*xl*yj);  // Xl Yj
						mkAsp(SS,2*(SS*(hno+l-1)+hno+l)-1,+GI*xl*yl);  // Xl Yl
						mkAsp(SS,2*(SS*(hno+l-1)+hno+m)-1,+GI*xl*ym);  // Xl Ym
		
						mkAsp(SS,2*(SS*(hno+m-1)+hno+j-1),+GI*xm*xj);  // Xm Xj
						mkAsp(SS,2*(SS*(hno+m-1)+hno+l-1),+GI*xm*xl);  // Xm Xl
						mkAsp(SS,2*(SS*(hno+m-1)+hno+m-1),+GI*xm*xm);  // Xm Xm
						mkAsp(SS,2*(SS*(hno+m-1)+hno+j)-1,+GI*xm*yj);  // Xm Yj
						mkAsp(SS,2*(SS*(hno+m-1)+hno+l)-1,+GI*xm*yl);  // Xm Yl
						mkAsp(SS,2*(SS*(hno+m-1)+hno+m)-1,+GI*xm*ym);  // Xm Ym
		
						mkAsp(SS,2*(SS*(hno+j)+hno+j-1)-SS,+GI*yj*xj);  // Yj Xj
						mkAsp(SS,2*(SS*(hno+j)+hno+l-1)-SS,+GI*yj*xl);  // Yj Xl
						mkAsp(SS,2*(SS*(hno+j)+hno+m-1)-SS,+GI*yj*xm);  // Yj Xm
						mkAsp(SS,2*(SS*(hno+j)+hno+j)-1-SS,+GI*yj*yj);  // Yj Yj
						mkAsp(SS,2*(SS*(hno+j)+hno+l)-1-SS,+GI*yj*yl);  // Yj Yl
						mkAsp(SS,2*(SS*(hno+j)+hno+m)-1-SS,+GI*yj*ym);  // Yj Ym
		
						mkAsp(SS,2*(SS*(hno+l)+hno+j-1)-SS,+GI*yl*xj);  // Yl Xj
						mkAsp(SS,2*(SS*(hno+l)+hno+l-1)-SS,+GI*yl*xl);  // Yl Xl
						mkAsp(SS,2*(SS*(hno+l)+hno+m-1)-SS,+GI*yl*xm);  // Yl Xm
						mkAsp(SS,2*(SS*(hno+l)+hno+j)-1-SS,+GI*yl*yj);  // Yl Yj
						mkAsp(SS,2*(SS*(hno+l)+hno+l)-1-SS,+GI*yl*yl);  // Yl Yl
						mkAsp(SS,2*(SS*(hno+l)+hno+m)-1-SS,+GI*yl*ym);  // Yl Ym
		
						mkAsp(SS,2*(SS*(hno+m)+hno+j-1)-SS,+GI*ym*xj);  // Ym Xj
						mkAsp(SS,2*(SS*(hno+m)+hno+l-1)-SS,+GI*ym*xl);  // Ym Xl
						mkAsp(SS,2*(SS*(hno+m)+hno+m-1)-SS,+GI*ym*xm);  // Ym Xm
						mkAsp(SS,2*(SS*(hno+m)+hno+j)-1-SS,+GI*ym*yj);  // Ym Yj
						mkAsp(SS,2*(SS*(hno+m)+hno+l)-1-SS,+GI*ym*yl);  // Ym Yl
						mkAsp(SS,2*(SS*(hno+m)+hno+m)-1-SS,+GI*ym*ym);  // Ym Ym
		
						bv[2*i-2]=bv[2*i-2]-GI*om*xi;  // Xi
						bv[2*i-1]=bv[2*i-1]-GI*om*yi;  // Yi
		
						bv[2*(hno+j-1)]=bv[2*(hno+j-1)]-GI*om*xj;  // Xj
						// bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-GI*om*xk;  // Xk
						bv[2*(hno+l-1)]=bv[2*(hno+l-1)]-GI*om*xl;  // Xl
						bv[2*(hno+m-1)]=bv[2*(hno+m-1)]-GI*om*xm;  // Xm
		
						bv[2*(hno+j)-1]=bv[2*(hno+j)-1]-GI*om*yj;  // Yj
						// bv[2*(hno+k)-1]=bv[2*(hno+k)-1]-GI*om*yk;  // Yk
						bv[2*(hno+l)-1]=bv[2*(hno+l)-1]-GI*om*yl;  // Yl
						bv[2*(hno+m)-1]=bv[2*(hno+m)-1]-GI*om*ym;  // Ym
		
						cconst=cconst+0.5L*GI*om*om;
					}
					else
					{
						// Case D
						if (((*h).st2)&&(((*(*h).st1).rod)!=((*(*h).st2).rod)))
						{
							// Case D: 2 sticks on different rods
							xi=-sg2/ll+sg1/lj;
							yi= cg2/ll-cg1/lj;
							xj=-sg1/lj;
							yj= cg1/lj;
							xl= sg2/ll;
							yl=-cg2/ll;
							om=(xj*dxj+yj*dyj+xl*dxl+yl*dyl);
							// Hinge-hinge types
							mkAsp(SS,2*(SS*(i-1)+i-1),+GI*xi*xi); // Xi Xi
							mkAsp(SS,2*(i*SS+i)-1-SS,+GI*yi*yi);  // Yi Yi
							mkAsp(SS,2*(SS*(i-1)+i)-1,+GI*xi*yi);  // Xi Yi
							mkAsp(SS,2*(i*SS+i-1)-SS,+GI*yi*xi);  // Yi Xi
			
							// Hinge-midpoint types
							mkAsp(SS,2*(SS*(i-1)+hno+j-1),+GI*xi*xj);  // Xi Xj
							mkAsp(SS,2*(SS*(i-1)+hno+l-1),+GI*xi*xl);  // Xi Xl
			
							mkAsp(SS,2*(i*SS+hno+j)-1-SS,+GI*yi*yj);  // Yi Yj
							mkAsp(SS,2*(i*SS+hno+l)-1-SS,+GI*yi*yl);  // Yi Yl
			
							mkAsp(SS,2*(SS*(i-1)+hno+j)-1,+GI*xi*yj);  // Xi Yj
							mkAsp(SS,2*(SS*(i-1)+hno+l)-1,+GI*xi*yl);  // Xi Yl
			
							mkAsp(SS,2*(i*SS+hno+j-1)-SS,+GI*yi*xj);  // Yi Xj
							mkAsp(SS,2*(i*SS+hno+l-1)-SS,+GI*yi*xl);  // Yi Xl
			
							// Midpoint-hinge types
							mkAsp(SS,2*(SS*(hno+j-1)+i-1),+GI*xj*xi);  // Xj Xi
							mkAsp(SS,2*(SS*(hno+l-1)+i-1),+GI*xl*xi);  // Xl Xi
			
							mkAsp(SS,2*(SS*(hno+j)+i)-1-SS,+GI*yj*yi);  // Yj Yi
							mkAsp(SS,2*(SS*(hno+l)+i)-1-SS,+GI*yl*yi);  // Yl Yi
			
							mkAsp(SS,2*(SS*(hno+j)+i-1)-SS,+GI*yj*xi);  // Yj Xi
							mkAsp(SS,2*(SS*(hno+l)+i-1)-SS,+GI*yl*xi);  // Yl Xi
			
							mkAsp(SS,2*(SS*(hno+j-1)+i)-1,+GI*xj*yi);  // Xj Yi
							mkAsp(SS,2*(SS*(hno+l-1)+i)-1,+GI*xl*yi);  // Xl Yi
			
							// Midpoint-midpoint types
							mkAsp(SS,2*(SS*(hno+j-1)+hno+j-1),+GI*xj*xj);  // Xj Xj
							mkAsp(SS,2*(SS*(hno+j-1)+hno+l-1),+GI*xj*xl);  // Xj Xl
							mkAsp(SS,2*(SS*(hno+j-1)+hno+j)-1,+GI*xj*yj);  // Xj Yj
							mkAsp(SS,2*(SS*(hno+j-1)+hno+l)-1,+GI*xj*yl);  // Xj Yl
			
							mkAsp(SS,2*(SS*(hno+l-1)+hno+j-1),+GI*xl*xj);  // Xl Xj
							mkAsp(SS,2*(SS*(hno+l-1)+hno+l-1),+GI*xl*xl);  // Xl Xl
							mkAsp(SS,2*(SS*(hno+l-1)+hno+j)-1,+GI*xl*yj);  // Xl Yj
							mkAsp(SS,2*(SS*(hno+l-1)+hno+l)-1,+GI*xl*yl);  // Xl Yl
			
							mkAsp(SS,2*(SS*(hno+j)+hno+j-1)-SS,+GI*yj*xj);  // Yj Xj
							mkAsp(SS,2*(SS*(hno+j)+hno+l-1)-SS,+GI*yj*xl);  // Yj Xl
							mkAsp(SS,2*(SS*(hno+j)+hno+j)-1-SS,+GI*yj*yj);  // Yj Yj
							mkAsp(SS,2*(SS*(hno+j)+hno+l)-1-SS,+GI*yj*yl);  // Yj Yl
			
							mkAsp(SS,2*(SS*(hno+l)+hno+j-1)-SS,+GI*yl*xj);  // Yl Xj
							mkAsp(SS,2*(SS*(hno+l)+hno+l-1)-SS,+GI*yl*xl);  // Yl Xl
							mkAsp(SS,2*(SS*(hno+l)+hno+j)-1-SS,+GI*yl*yj);  // Yl Yj
							mkAsp(SS,2*(SS*(hno+l)+hno+l)-1-SS,+GI*yl*yl);  // Yl Yl
			
							bv[2*i-2]=bv[2*i-2]-GI*om*xi;  // Xi
							bv[2*i-1]=bv[2*i-1]-GI*om*yi;  // Yi
			
							bv[2*(hno+j-1)]=bv[2*(hno+j-1)]-GI*om*xj;  // Xj
							// bv[2*(hno+k-1)]=bv[2*(hno+k-1)]-GI*om*xk;  // Xk
							bv[2*(hno+l-1)]=bv[2*(hno+l-1)]-GI*om*xl;  // Xl
							// bv[2*(hno+m-1)]=bv[2*(hno+m-1)]-GI*om*xm;  // Xm
			
							bv[2*(hno+j)-1]=bv[2*(hno+j)-1]-GI*om*yj;  // Yj
							// bv[2*(hno+k)-1]=bv[2*(hno+k)-1]-GI*om*yk;  // Yk
							bv[2*(hno+l)-1]=bv[2*(hno+l)-1]-GI*om*yl;  // Yl
							// bv[2*(hno+m)-1]=bv[2*(hno+m)-1]-GI*om*ym;  // Ym
			
							cconst=cconst+0.5L*GI*om*om;
						}
						else
						{
							// This case is not pissble
							printf("IMPOSSIBLE CASE IN SHEAR!!!\n");
						}
					}
				}
			}
		}
		h=(*h).next;
	}
	return 0;
}

// **************************

int conjgrad(int hno, int mpno)
// Conjugate gradient method to find the solution to Ax=b, the solution x
// gives minimizer of energy H(x)=xAx/2-bx+c.
{
	int SS,i,iter;
	long double alpha,beta,nev,r2;
	long double *r,*p,*Ap,*xvec;
	struct hinges *hh;
	struct midpoints *mm;
	SS=2*(hno+mpno);
	xvec=(long double *)calloc(SS,sizeof(long double));
	r=(long double *)calloc(SS,sizeof(long double));
	Ap=(long double *)calloc(SS,sizeof(long double));
	p=(long double *)calloc(SS,sizeof(long double));
	// Initial guess for solution is current configuration: load to vector xvec
	i=-1;
	hh=hingestart;
	while (hh)
	{
		i++;
		xvec[i]=(*hh).x;
		i++;
		xvec[i]=(*hh).y;
		hh=(*hh).next;
	}
	mm=midpointstart;
	while (mm)
	{
		i++;
		xvec[i]=(*mm).x;
		i++;
		xvec[i]=(*mm).y;
		mm=(*mm).next;
	}
	// Now xvec contains the actual configuration. Compute residual:
	for (i=0;i<SS;i++)
	{
		r[i]=bv[i];
	}
	for (i=1;i<=kspr;i++)
	{
		r[Aspri[i]-1]=r[Aspri[i]-1]-Aspr[i]*xvec[Asprj[i]-1];
	}
	r2=0.0L;
	for (i=0;i<SS;i++)
	{
		p[i]=r[i];
		r2=r2+r[i]*r[i];
	}

// Cycle until we reach solution
	iter=0;
//	while ((r2>SMALL)&&(iter<SS))  //   ÉRDEKES, HOGY NEM KONVERGAL SS LÉPÉSBEN. MIÉRT? KÉNYSZEREK MIATT? EMIATT INKÁBB A KÖVETKEZŐ SOR LETT:
	while ((r2>SMALL))
	{
		iter++;
		nev=0.0L;
		for (i=0;i<SS;i++)
		{
			Ap[i]=0.0L;
		}
		for (i=1;i<=kspr;i++)
		{
			Ap[Aspri[i]-1]=Ap[Aspri[i]-1]+p[Asprj[i]-1]*Aspr[i];
		}
		nev=0.0L;
		for (i=0;i<SS;i++)
		{
			nev=nev+Ap[i]*p[i];
		}
		alpha=r2/nev;
		beta=0.0L;
		for (i=0;i<SS;i++)
		{
			r[i]=r[i]-alpha*Ap[i];
			beta=beta+r[i]*r[i];
		}
		beta=beta/r2;
		r2=0.0L;
		for (i=0;i<SS;i++)
		{
			xvec[i]=xvec[i]+alpha*p[i];
			p[i]=r[i]+beta*p[i];
			r2=r2+r[i]*r[i];
		}
	}
//printf("SS=%d, ITER= %d, r2= %20.16Lf\n",SS,iter,r2);
	// Solution is in xvec, write it back to network data
	i=-1;
	hh=hingestart;
	while (hh)
	{
		i++;
		(*hh).x=xvec[i];
		i++;
		(*hh).y=xvec[i];
		hh=(*hh).next;
	}
	mm=midpointstart;
	while (mm)
	{
		i++;
		(*mm).x=xvec[i];
		i++;
		(*mm).y=xvec[i];
		mm=(*mm).next;
	}
	return 0;
}

// ***********************
// Main program
// ***********************

int main()
{
	int err,hno,mpno,cyc;
	char fname[30];
	long double energy1,energy2,energy3,energy;
	time_t rawtime;
	struct tm * timeinfo;
	FILE * statusfile;

	// Set up random seed
	srand48(SEED);
	// Initiate files
	err=setupfiles();

	sprintf(fname,"status%04d.dat",FILENUMBER);
	statusfile = fopen(fname,"w");
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        fprintf (statusfile, "%d - program started at: %s", FILENUMBER, asctime (timeinfo) );
        printf ( "%d - program started at: %s", FILENUMBER, asctime (timeinfo) );
	fclose(statusfile);

	// Cycle along all networks with same data
	for (cyc=0;cyc<CASENO;cyc++)
	{
		// Setup rods
		err=setuprods();
		// Print rods
//		err=printrods();
		// Look for intersections between rods, hno is number of hinges found
		hno=findhinges();
		// Find all sticks (and fill in sticks, midpoints, and from hinges: st1,2,3,4 n1,2,3,4 mp1,2,3,4)
		// mpno is number of midpoints (and hence sticks) found
		mpno=findsticks();
		// Delete dangling ends, that is, those hinges which has no attached stick. hno is final number of hinges.
		hno=deletehinges();
		// Print base data for the network
		err=networkdata();

		// Prints hinges, sticks, midpoints and energy before loading
//		printf("ENERGIES BEFORE LOADING\n");
//		err=printhinges();
//		err=printmidpoints();
//		err=printsticks(0.0L,0.0L);
//		energy1=stretchenergy(0.0L,0.0L);
//		energy2=bendenergy(0.0L,0.0L);
//		energy3=shearenergy(0.0L,0.0L);
//		printf("ENERGIES BEFORE LOAD\n%20.16Lf  %20.16Lf  %20.16Lf  %20.16Lf\n",energy1,energy2,energy3,energy1+energy2+energy3);
//		printf("=========================\n\n");

		// Prints hinges, sticks, midpoints and energy after the loading, but before relaxation
//		printf("ENERGIES AFTER LOADING, BEFORE BALANCE\n");
//		err=printhinges();
//		err=printmidpoints();
//		err=printsticks(SHEARDISPL,STRETCHDISPL);
//		energy1=stretchenergy(SHEARDISPL,STRETCHDISPL);
//		energy2=bendenergy(SHEARDISPL,STRETCHDISPL);
//		energy3=shearenergy(SHEARDISPL,STRETCHDISPL);
//		printf("ENERGIES AFTER UNBALANCED LOAD\n%20.16Lf  %20.16Lf  %20.16Lf  %20.16Lf\n",energy1,energy2,energy3,energy1+energy2+energy3);
	
		// Look for solution (equilibrium) by minimising energy xAx/2-bx+c
		// First define vector b
		bv=(long double *)calloc(2*(hno+mpno),sizeof(long double));
		// Load matrix A and vector b for conjugate gradient method
		err=conjgradvars(hno,mpno,SHEARDISPL,STRETCHDISPL);
	
		// Test: print energy from matrix: xAx/2-bx+c
//		energy=energyfroma(hno,mpno);
//		printf("Energy from A: %20.16Lf\n",energy);
//		printf("DIFFERENCE:    %20.16Lf\n",energy-energy1-energy2-energy3);
//		printf("=========================\n\n");
	
		// Sequence of lines and rows in A and elemenst in b are the same as the
		// sequence of elements in chain 'hinges' then elements in chain 'midpoints'
		// Look for solution
		err=conjgrad(hno,mpno);

		// Prints hinges, sticks, midpoints and energy after relaxation (in equilibrium
//		printf("ENERGIES AFTER BALANCING:\n");
//		err=printhinges();
//		err=printmidpoints();
//		err=printsticks(SHEARDISPL,STRETCHDISPL);
		energy1=stretchenergy(SHEARDISPL,STRETCHDISPL);
		energy2=bendenergy(SHEARDISPL,STRETCHDISPL);
		energy3=shearenergy(SHEARDISPL,STRETCHDISPL);
//		printf("ENERGIES IN BALANCE\n%20.16Lf + %20.16Lf + %20.16Lf = %20.16Lf\n",energy1,energy2,energy3,energy1+energy2+energy3);
		energy=energyfroma(hno,mpno);
//		printf("Energy from A: %20.16Lf\n",energy);
//		printf("DIFFERENCE:    %20.16Lf\n",energy-energy1-energy2-energy3);
//		printf("=========================\n\n");
	
		// Write final energy in data file
		fprintf(data,"%20.16Lf %20.16Lf %20.16Lf %20.16Lf %20.16Lf %20.16Lf\n",energy1,energy2,energy3,energy1+energy2+energy3,energy,(energy1+energy2+energy3-energy)/energy);
		
		// Clean up old network
		err=cleanuprods();
		err=cleanuphinges();
		err=cleanupmidpoints();
		err=cleanupsticks();
		free(bv);
		free(Aspr);
		free(Aspri);
		free(Asprj);
		statusfile = fopen(fname,"a");
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf (statusfile, "%d - cycle %d ended at: %s", FILENUMBER, cyc, asctime (timeinfo) );
		printf ( "%d - cycle %d ended at: %s", FILENUMBER, cyc, asctime (timeinfo) );
		fclose(statusfile);
	}
	// End of program
	fclose(data);
	fclose(roddata);
	fclose(hingedata);
	fclose(midpointdata);
	fclose(stickdata);
	err=0;
	return err;
}

