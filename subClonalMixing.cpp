
# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <getopt.h>



using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 1e5;

// Next-nearest-neighbour neighbourhood -> 8 neighbours 
const int NEIGHBOURHOOD = 8;

double s, tmut;
int q, seed;

double radius_double, t, r_birth_WT, r_birth_mutant, r_birth_WT_normalised, r_birth_mutant_normalised, rand_double;
double optimal_direction_i,  optimal_direction_j, optimal_vector_norm, vector_norm, rescaled_min_length, scalar_prod;
int radius, Ntot, Nwt, iter, x, y, cell_x, cell_y, dir, queue, ind, length, coordX, coordY, previous_link_direction;
int arising_time, x_b, y_b, direction, chosen_direction, min_length, num_mins, chosen_min;



// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell 
int chainX[(int)_maxsize];
int chainY[(int)_maxsize];
int directions[NEIGHBOURHOOD];



bool mutated = false;
bool extended_chain = false;
bool quiet = true;
bool chainstuck = false;
bool WT_BIRTH , MUTANT_BIRTH;



/*******************************************************************************/



// Define poisson distributions
default_random_engine generator;



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int dvr;

	// Constructor for Cell object
	Cell(){}

	// Set() and get() methods
	
	void setDVR(int n)
	{
		this->dvr = n;
	}
	

	void setGAs(int d)
	{
		this->dvr = d;
	}

};







// Surface growth with division rate proportional to number of empty neighbours
void surface_division(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nwt , int *x_b , int *y_b , int radius )
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	if (tumour[cell_x + x][cell_y + y].dvr == -1)		// Check if neighbour is empty
	{


		// Create daughter cell
		tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr);

		*Ntot += 1;

		if (tumour[cell_x + x][cell_y + y].dvr == 0)
		{
			*Nwt += 1;
		}


		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
	}

}









// Define cell division in volumetric growth model as described in "Waclaw, B., et al., Nature 525, 261–264 (2015)" with addition of "Newton 2" cell displacement
void volumetric_division(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nwt , int *x_b , int *y_b , int radius )
{





	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));



	// Find shortest path to an empty lattice point in the tumour
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;

	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;

	while(1)
	{

		if (tumour[chainX[queue]][chainY[queue]].dvr == -1) break;


		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

					if ( (i == 0) && (j == 0) )
					{
						//++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;


					if (directions[direction] == _maxsize)
					{
						//cout << "Skipping (dx,dy,dz) = (" << i << " , " << j << " , " << k <<  ") --> direction=" << directions[direction] << endl;
						++direction;
						continue;
					}
 
					//cout << "Searching direction (dx,dy,dz) = (" << i << " , " << j << " , " << k <<  ")" << endl;

					length = 1;
					while(1)
					{
						coordX = chainX[queue] + length*i;
						coordY = chainY[queue] + length*j;

						//cout << "The status of (dx,dy,dz) = (" << coordX << " , " << coordY  << " , " << coordZ << "   is:  " << tumour[coordX][coordY][coordZ].cellID << endl;

						if (tumour[coordX][coordY].dvr == -1)
						{
							directions[direction] = length;
							break;
						}
						else ++length;
					}

					++direction;
			}
		}

		//cout << "Distance to empty cell in each direction: ";
		//for (int i = 0; i < NEIGHBOURHOOD; ++i)
		//{
		//	cout << directions[i] << " ";
		//}
		//cout << endl;

		extended_chain = false;
		while(extended_chain == false)
		{


			// Find which entry in directions list is smallest
			min_length = *Ntot;
			num_mins = 0;
			chainstuck = true;
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
				if (directions[i] <= _maxsize) chainstuck = false;
			}

			if (chainstuck == true) 
			{
				//cout << "Chain stuck" << endl
				return;
			}



			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}

			//cout << "Number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;


			
			if (queue >= 1)
			{
				// Check in which direction the previous link in the chain is at
				direction = 0;
				previous_link_direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if ( (chainX[queue]+i == chainX[queue-1]) && (chainY[queue]+j == chainY[queue-1]) )
						{
							previous_link_direction = direction;


							//cout << "Direction of previous link in chain: (" << chainX[queue-1] << "," << chainY[queue-1] << ") -> (" << chainX[queue] << "," << chainY[queue] << ")  ====  " << previous_link_direction << ", opposite direction -> " << 7-previous_link_direction << endl;
							break;
						}

						++direction;

					}
				}

				// Find direction the cell would travel in if unrestricted by considering the which direction it is being pushed from, i.e. where it's neighbour is relative to itself (Newton's 2nd law)
				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (direction == NEIGHBOURHOOD - 1 - previous_link_direction)	// Opposite direction
						{
							optimal_direction_i = (double)i;
							optimal_direction_j = (double)j;
							optimal_vector_norm = pow(((optimal_direction_i*optimal_direction_i) + (optimal_direction_j*optimal_direction_j)) , 0.5);

							//cout << "(i,j)=(" << i << "," << j << ") -> normalising factor = " << optimal_vector_norm << endl;

							optimal_direction_i /= optimal_vector_norm;
							optimal_direction_j /= optimal_vector_norm;

							//cout << "Optimal direction = " << optimal_direction_i << " " << optimal_direction_j << endl;
						}

						++direction;

					}
				}


				// Re-scale distances vector according to relative direction to 'forward' chain direction
				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						vector_norm = pow(((i*i) + (j*j)) , 0.5);

						// Scalar product with unit vector pointing in 'optimal direction'
						scalar_prod = (optimal_direction_i*i/vector_norm) + (optimal_direction_j*j/vector_norm);

						//cout << "Candidate direction = " << i/vector_norm << " " << j/vector_norm << endl;

						// Rescale to within range [0,1]
						scalar_prod = (scalar_prod + 1.0)/2.0;
						//cout << "Mutliplying " << directions[direction] << " by " << 1.0 - scalar_prod << endl;
						directions[direction] *= 1.0 - scalar_prod;
						//directions[direction] /= scalar_prod;

						++direction;

					}
				}



				// Find new minimum after rescaling 
				rescaled_min_length = (double)*Ntot;
				num_mins = 0;
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] < min_length) min_length = directions[i];
				}

				// Then count number of directions which are minimum
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) ++num_mins;
				}

				//cout << "After rescaling: distance to empty cell in each direction: ";
				//for (int i = 0; i < NEIGHBOURHOOD; ++i)
				//{
				//	cout << directions[i] << " ";
				//}
				//cout << endl;

				//cout << "After rescaling, new number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;

			}
			

			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{


				ind = 0;
				
				while(1)
				{
					rand_double = drand48();
					//if (directions[ind%NEIGHBOURHOOD] == min_length) 
					//{
					//	cout << "Checking directions[" << ind%NEIGHBOURHOOD << "]" << endl;
					//	cout << "Checking if " << ran << "<" << 1.0/(float)num_mins << endl;
					//}
					if ( (directions[ind%NEIGHBOURHOOD] == min_length) && (rand_double < (1.0/(double)num_mins)))
					{
						chosen_direction = ind%NEIGHBOURHOOD;
						break;
					}



					++ind;
				}
			}

			else 	// Otherwise select the only minimal direction
			{
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) chosen_direction = i;
				}
			}


			//cout << "Chosen direction is: " << chosen_direction << endl;


			// Before adding the next cell in the chosen direction to the chain, make sure chain is self-avoiding
			// First, find the coordinates of potential new cell 
			direction = 0;
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
					if ( (i == 0) && (j == 0) ) 
					{
						continue;
					}
					//cout << direction << endl;

					if (direction == chosen_direction)
					{
						coordX = chainX[queue] + i;
						coordY = chainY[queue] + j;

						//cout << "Previous chain coord (x,y,z) = " << chainX[queue] << " , " << chainY[queue] << ")" << endl;
						//cout << "Next link at coord (x,y,z) = " << coordX << " , " << coordY << " , " << ")" << endl;

					}

					++direction;

				}
			}

			// Second, check these coordinates are not already in the chain
			extended_chain = true;
			for (int i = 0; i < (queue+1); ++i)
			{

				//cout << "Checking queue entry " << i << "/" << queue << endl;
				//if ( ((chainX[i] == coordX) && (chainY[i] == coordY) && (chainZ[i] == coordZ)) && ((coordX == cell_x) && (coordY == cell_y) && (coordZ == cell_z)) )
				if ( ((chainX[i] == coordX) && (chainY[i] == coordY)) || ((coordX == cell_x) && (coordY == cell_y)) )
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					//cout << "Potential new link in chain is already taken! ****** " << extended_chain << endl;
					directions[chosen_direction] += (int)_maxsize;
					extended_chain = false;
				}
			}

			if (extended_chain == true)
			{
				//cout << "Potential new link in chain is not taken! ******" << endl;
			}


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY;
		

		++queue;


		// Put a limit on the number of cells that any particular cell can push
		if (queue > q)
		{
			//cout << "Queue exceeded" << endl;
			return;
		}

	}




	// If cell can divide (i.e. not in quiescent state) then cell dies with probability p_death = 1/3
	if ((drand48() < 1.0/3.0) && !(*Ntot == 1))
	{

		//cout << "cell died during division" << endl;


		*Ntot -= 1;
		if (tumour[cell_x][cell_y].dvr == 0)
		{
			*Nwt -= 1;
		}


		tumour[cell_x][cell_y].dvr = -1;

		return;
	}





	// Once the chain has been constructed, move all cells along one place
	if (queue > 0)
	{

		for (int i = 0; i < queue; ++i)
		{
			//cout << "chain(x" << i << ",y" << i << ") = (" << chainX[i] << " , " << chainY[i] << ") -> " << tumour[chainX[i]][chainY[i]].dvr << " ||| Pushing cell at (x,y)=(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") to (x,y)=(" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;
			//tumour[chainX[queue-i]][chainY[queue-i]][chainZ[queue-i]].setGAs(tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].dvr , tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].res , tumour[chainX[queue-i-1]][chainY[queue-i-1]][chainZ[queue-i-1]].pgr);
			tumour[chainX[queue-i]][chainY[queue-i]].dvr = tumour[chainX[queue-i-1]][chainY[queue-i-1]].dvr;


			// Update bounds on tumour size
			if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
			if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);

			
		}
	}

	else 			// Even if queue=0, check that newly created cell increases any bounds
	{

		// Update bounds on tumour size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);

	}


	
	// Create daughter cell
	tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr);


	*Ntot += 1;


	if (tumour[cell_x + x][cell_y + y].dvr == 0)
	{
		*Nwt += 1;
	}


	//cout << "cell divided" << endl;

			
}










/*******************************************************************************/










int main(int argc, char** argv)
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();



	// Reset time and tumour size variables
	t = 0.0;
	Ntot = 0;
	Nwt = 0;




	// Parse command line arguments
	int c;

	while ((c = getopt (argc, argv, ":v:x:s:t:q:")) != -1)
	switch (c)
	{
		// Verbose flag
		case 'v':
			quiet = false;
			break;

		case 'x':
			seed = atoi(optarg);		
			break;

		// Set selective advantage of mutant cells
		case 's':
			s = atof(optarg);		
			break;

		// Set time of emergence of mutant clone, as a percentage of max size
		case 't':
			tmut = atof(optarg);		
			break;

		// Set pushing limit for division algorithm
		case 'q':
			q = atoi(optarg);		
			break;

		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		return 1;
		default:
		abort ();
	}



	// Seed random number generator
	srand48(seed);

	
	// Simple checks in mutation time parameter
	if ((tmut < 0.0) || (tmut > 1.0))
	{
		cout << "Mutation time parameter, t, must be given as a fraction of final system size i.e. 0 <= t <= 1" << endl;
		exit(0);
	}
	arising_time = (int)(tmut * (double)(_maxsize));








	//================== Initialise tumour ====================//

	// Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
	radius_double = pow ( (_maxsize/M_PI) , (1.0/2.0) );
	radius = (int)(4*radius_double);

	if (!quiet) cout << " " << endl;

	Cell ** tumour = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		tumour[i] = new Cell[2*radius];

		for (int j = 0; j < (2*radius); j++)
		{

			tumour[i][j].setDVR(-1);		// Cells with dvr = -1 represent empty lattice points, dvr = 0 and dvr = 1 represent WT and mutant cells respectively
		
		}
		if (!quiet) printf(" Initialising tumour... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}
	if (!quiet) printf(" Initialising tumour... Done.\r");
	if (!quiet) cout << " " << endl;
		


	// Seed first tumour cell/gland at (x,y) = (radius , radius)
	tumour[radius][radius].setDVR(0);


	Ntot += 1;
	Nwt += 1;










	//================== Simulate tumour growth ==================//

	iter = 0;
	x = 0;
	y = 0;
	x_b = 10;
	y_b = 10;


	do
	{
		
		++iter;

		//cout << "N = " << Ntot << endl;


		// Choose reaction (PT and non-PT cell division considered different reactions)
		r_birth_WT = 1.0;
		r_birth_mutant = r_birth_WT + s;



		// Multiply reaction rates by the relevant number of cells
		r_birth_WT *= Nwt;
		r_birth_mutant *= Ntot - Nwt;



		// Compute normalised reaction rates
		r_birth_WT_normalised = r_birth_WT/(r_birth_WT + r_birth_mutant);
		r_birth_mutant_normalised = r_birth_mutant/(r_birth_WT + r_birth_mutant);



		if (Ntot > 0.1*_maxsize)
		{
			x_b = radius;
			y_b = radius;
		}





		// Create mutant cell at specified time
		if ((Ntot == arising_time) && !(mutated))
		{
			mutated = true;

			//cout << "\nMutation appeared...\n" << endl;

			// Randomly select cell to mutate
			do
			{
				cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
				cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
			}
			while (tumour[cell_x][cell_y].dvr == -1);

			tumour[cell_x][cell_y].dvr = 1;
			Nwt -= 1;
		}






		// If mutant clone dies out, exit
		if ((Ntot > arising_time) && (Ntot-Nwt == 0))
		{
			cout << "0" << endl;
			exit(0);
		}






		// Gillespie rates
		WT_BIRTH = false;
		MUTANT_BIRTH = false;

		rand_double = drand48();
		if (rand_double <= r_birth_WT_normalised)
		{
			WT_BIRTH = true;
		}
		else if (rand_double <= r_birth_WT_normalised + r_birth_mutant_normalised)
		{
			MUTANT_BIRTH = true;
		}
		else
		{
			cout << "Problem with Gillespie rates..." << endl;
			exit(0);
		}







		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;

		do
		{

			cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
		}
		while (((WT_BIRTH) && (tumour[cell_x][cell_y].dvr != 0)) || ((MUTANT_BIRTH) && (tumour[cell_x][cell_y].dvr != 1)));


		//cout << "Chosen cell is (x,y) = (" << cell_x << " , " << cell_y << ")" << endl;


		// Cell divides (or may die)
		volumetric_division(tumour , cell_x , cell_y , &Ntot , &Nwt , &x_b , &y_b , radius);








		if (iter%10000 == 0)
		{
			if (!quiet) cout << "Iter #" << iter << ": Ntot=" << Ntot << " , Nwt=" << Nwt << " , Nmut=" << Ntot-Nwt << " ||| Bounds: (x) " << radius - x_b << "--" << radius + x_b << ", (y) " << radius - y_b << "--" << radius + y_b << endl;
		}
			


	} while (Ntot < _maxsize);		// Exit once system has reached total size of _maxsize

	if (!quiet) cout << " " << endl;




	// Systems with a final mutant subclone size of <200 cells are deemed to be unsuccessful
	if (Ntot-Nwt < 200)
	{
		cout << "0" << endl;
		exit(0);
	}







	//================== Open data files ==================//

	stringstream f;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_s=" << s << "_tmut=" << tmut << "_q=" << q << "/seed=" << seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./2D_DATA/maxSize=" << _maxsize << "_s=" << s << "_tmut=" << tmut << "_q=" << q << "/seed=" << seed;
		system(f.str().c_str());
	}

	ofstream tumour_file;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_s=" << s << "_tmut=" << tmut << "_q=" << q << "/seed=" << seed << "/cells.csv";
	tumour_file.open(f.str().c_str());


	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "Created output files..." << endl;





	// Write tumour data to file
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			if ( tumour[i][j].dvr != -1 ) tumour_file << i << "," << j << "," << tumour[i][j].dvr << endl;
		}
	}



	tumour_file.close();


	cout << "1" << endl; 

	return 0;
}

















