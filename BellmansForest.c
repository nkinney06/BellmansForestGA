#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

typedef struct {
    float x;
    float y;
} Cartesian;

typedef struct {
    float r;
    float theta;
} Polar;

// pair of point lists for the escape path
Cartesian randomWalkCartesian[100];
Polar randomWalkPolar[100];
Polar population[1000][100];
Polar newPopulation[1000][100];

// points for the forest and grid
Cartesian points[10000];
Cartesian forest[10000];

// default parameters of the system
int degreesOfFreedom = 2;
int totalForestSides = 3;
float forestGridSize = 5.;
int populationSize = 100;
int generations = 100;
int children = 40;
int survivors = 40;

// additional globals
int numberTestPoints = 0;
int numberTestAngles = 100;
int worsteAngle;
int worstePoint;
float fitness[1000];
int minFitnessIndex = 0;

/////////////////////////////////////////////////////////////////////////////////////
// initialization functions
/////////////////////////////////////////////////////////////////////////////////////
void makeWalk() {
	for ( int i = 0 ; i < degreesOfFreedom ; i++ ) {
		randomWalkPolar[i].theta = ((float)rand() / RAND_MAX) * 2.0f * M_PI;
		randomWalkPolar[i].r = 100.0 * rand() / (float) RAND_MAX;
	}
	randomWalkPolar[0].theta = 0;
	randomWalkPolar[degreesOfFreedom-1].r = 300;
}

void clearCartesian() {
	for ( int i = 0 ; i < 100 ; i++ ) {
		randomWalkCartesian[i].x = 0;
		randomWalkCartesian[i].y = 0;
	}
}

void makeForest() {
	float t = 0;
	float r = 0;
	for ( int i = 0 ; i < totalForestSides ; i++ ) {		
		t = i * 2.0f * M_PI / totalForestSides + ((float)rand() / RAND_MAX) * 2.0f * M_PI / totalForestSides;
		r = 100.0 * rand() / (float) RAND_MAX;
		forest[i].x = cos(t) * r;
		forest[i].y = sin(t) * r;
	}
}

// determine the orientation of three points (p, q, r)
int orientation(Cartesian p, Cartesian q, Cartesian r) {
    float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0.0f) return 0; // Collinear
    return (val > 0.0f) ? 1 : 2; // Clockwise or counterclockwise
}

// Function to test if a point is inside a polygon using the winding number algorithm
int isInsideForest(Cartesian testPoint) {
    int windingNumber = 0;
    for (int i = 0; i < totalForestSides; ++i) {
        Cartesian p1 = forest[i];
        Cartesian p2 = forest[(i + 1) % totalForestSides]; // Wrap around for the last edge
        if (p1.y <= testPoint.y) {
            if (p2.y > testPoint.y && orientation(p1, p2, testPoint) == 1) {
                windingNumber++;
            }
        } else {
            if (p2.y <= testPoint.y && orientation(p1, p2, testPoint) == 2) {
                windingNumber--;
            }
        }
    }
    return windingNumber != 0;
}

void gridForest() {
	Cartesian p;
	numberTestPoints = 0;
	for ( p.x = -100 ; p.x < 100 ; p.x = p.x + forestGridSize )
		for ( p.y = -100 ; p.y < 100 ; p.y = p.y + forestGridSize )
			if ( isInsideForest(p) )
				if (numberTestPoints < 10000)
					points[numberTestPoints++] = p;
}

void readPointsFromFile(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL)
        return;
    totalForestSides = 0;  // Number of points read
    while (fscanf(file, "%f\t%f", &forest[totalForestSides].x, &forest[totalForestSides].y) == 2)
        totalForestSides++;
    fclose(file);
}

/////////////////////////////////////////////////////////////////////////////////////
// conversion functions
/////////////////////////////////////////////////////////////////////////////////////
void cartesianToPolar() {
    Polar result;
	for (int i = 2; i < degreesOfFreedom; i++) {
		result.r = sqrt((randomWalkCartesian[i].x - randomWalkCartesian[i-1].x) * (randomWalkCartesian[i].x - randomWalkCartesian[i-1].x) + 
		                (randomWalkCartesian[i].y - randomWalkCartesian[i-1].y) * (randomWalkCartesian[i].y - randomWalkCartesian[i-1].y));
		result.theta = atan2(randomWalkCartesian[i].y - randomWalkCartesian[i-1].y, randomWalkCartesian[i].x - randomWalkCartesian[i-1].x);
		randomWalkPolar[i-1] = result;
	}
	randomWalkPolar[0].r = sqrt( pow( randomWalkCartesian[1].x-randomWalkCartesian[0].x,2 ) 
	                           + pow( randomWalkCartesian[1].y-randomWalkCartesian[0].y,2 ) );
	randomWalkPolar[0].theta = atan(randomWalkCartesian[1].y/randomWalkCartesian[1].x);  
}

void polarToCartesian() {
	clearCartesian();
    for (int i = 0; i < degreesOfFreedom; i++) {
        randomWalkCartesian[i+1].x = randomWalkCartesian[i].x + randomWalkPolar[i].r * cos(randomWalkPolar[i].theta);
        randomWalkCartesian[i+1].y = randomWalkCartesian[i].y + randomWalkPolar[i].r * sin(randomWalkPolar[i].theta);
    }
}

void translateCartesian(float x, float y) {
	for (int i = 0; i < degreesOfFreedom + 1; i++) {
		randomWalkCartesian[i].x = randomWalkCartesian[i].x + x;
		randomWalkCartesian[i].y = randomWalkCartesian[i].y + y;
	}
}

void rotateCartesian(float angle) {
	for (int i = 1; i < degreesOfFreedom + 1; i++) {
		float x = randomWalkCartesian[i].x;
		float y = randomWalkCartesian[i].y;
		randomWalkCartesian[i].x = x * cos(angle) - y * sin(angle);
		randomWalkCartesian[i].y = x * sin(angle) + y * cos(angle);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
// geometry functions
/////////////////////////////////////////////////////////////////////////////////////

// https://cboard.cprogramming.com/c-programming/154196-check-if-two-line-segments-intersect.html
// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in intersection.
char intersectSegments(Cartesian p0, Cartesian p1, Cartesian p2, Cartesian p3, Cartesian *intersection) {
    float s02_x, s02_y, s10_x, s10_y, s32_x, s32_y, s_numer, t_numer, denom, t;
    s10_x = p1.x - p0.x;
    s10_y = p1.y - p0.y;
    s02_x = p0.x - p2.x;
    s02_y = p0.y - p2.y;
 
    s_numer = s10_x * s02_y - s10_y * s02_x;
    if (s_numer < 0)
        return 0; // No collision
 
    s32_x = p3.x - p2.x;
    s32_y = p3.y - p2.y;
    t_numer = s32_x * s02_y - s32_y * s02_x;
    if (t_numer < 0)
        return 0; // No collision
 
    denom = s10_x * s32_y - s32_x * s10_y;
    if (s_numer > denom || t_numer > denom)
        return 0; // No collision
 
    // Collision detected
    t = t_numer / denom;
    intersection->x = p0.x + (t * s10_x);
	intersection->y = p0.y + (t * s10_y);
    return 1;
}

int detectEscape(Cartesian P, Cartesian Q, Cartesian* intersection) {
	float distance;
	float minDistance = 1000.0;
	int escapeSegment = -1;
	
	for (int i = 0; i < totalForestSides; i++) { // look for min distance to edge of forest
		if (intersectSegments(P,Q,forest[i],forest[(i+1)%totalForestSides],intersection)) {
			distance = sqrt( pow( P.x - intersection->x,2 ) + pow( P.y - intersection->y,2 ) );
			if ( distance < minDistance ) {
				escapeSegment = i;
				minDistance = distance;
			}
		}
	}
	
	if ( escapeSegment >= 0 ) { // if collision is detected store the min intersection point and return 1
		intersectSegments(P,Q,forest[escapeSegment],forest[(escapeSegment+1)%totalForestSides],intersection);
		return 1;
	}
	
	return 0;
}

float worsteCaseScenario(int *maxAngle, int *maxPoint) {
	Cartesian intersection;
	float distance,worseCase;
	int testPoint,segment,testAngle;

	worseCase = 0; // worste case scenario for this strat
	for ( testAngle = 0 ; testAngle < numberTestAngles; testAngle++ ) {
		for (testPoint = 0; testPoint < numberTestPoints; testPoint++) {
			translateCartesian(points[testPoint].x,points[testPoint].y);
			distance = 0; // walk along the escape path and add up the distance
			for (segment = 0; segment < degreesOfFreedom; segment++) {
				if ( detectEscape(randomWalkCartesian[segment],randomWalkCartesian[segment+1],&intersection) ) {
					distance += sqrt( pow( randomWalkCartesian[segment].x-intersection.x,2 ) + pow( randomWalkCartesian[segment].y-intersection.y,2 ) );
					if ( distance > worseCase ) {
						worseCase = distance;
						*maxAngle = testAngle;
						*maxPoint = testPoint;
					}
					break;
				} else { // still in the forest
					distance += randomWalkPolar[segment].r;
				}
			}
			translateCartesian(-1*points[testPoint].x,-1*points[testPoint].y);
		}
		rotateCartesian(2 * M_PI / numberTestAngles);
	}
	return worseCase;
}


/////////////////////////////////////////////////////////////////////////////////////
// output functions
/////////////////////////////////////////////////////////////////////////////////////
void printCartesian() {
	for (int i = 0; i < degreesOfFreedom + 1; i++)
		printf("%.2f\t%.2f\tstrats\n", randomWalkCartesian[i].x, randomWalkCartesian[i].y);
}

void printPolar() {
	printf("\n\npolar\n");
	for (int i = 0; i < degreesOfFreedom; i++)
		printf("Step %d: (%.2f, %.2f)\n", i, randomWalkPolar[i].r, randomWalkPolar[i].theta);
}

void printForest() {
	for (int i = 0; i < totalForestSides; i++)
		printf("%.2f\t%.2f\tforest\n", forest[i].x, forest[i].y);
}

void printGrid() {
	for (int i = 0; i < numberTestPoints; i++)
		printf("%.2f\t%.2f\tpoints\n", points[i].x, points[i].y);
}

void printRotations() {
	for ( int testAngle = 0 ; testAngle < numberTestAngles; testAngle++ ) {
		rotateCartesian(2 * M_PI / numberTestAngles);
		for (int i = 0; i < degreesOfFreedom + 1; i++)
			printf("%.2f\t%.2f\tcircle\n",randomWalkCartesian[i].x,randomWalkCartesian[i].y);
	}
}

float printSolution(int worsteAngle,int worstePoint) {
	Cartesian intersection;
	float distance = 0;
	rotateCartesian(worsteAngle * 2 * M_PI / numberTestAngles);
	translateCartesian(points[worstePoint].x,points[worstePoint].y);
	printf("%.2f\t%.2f\tsolved\n",randomWalkCartesian[0].x,randomWalkCartesian[0].y);
	for (int segment = 0; segment < degreesOfFreedom; segment++) {
		if ( detectEscape(randomWalkCartesian[segment],randomWalkCartesian[segment+1],&intersection) ) {
			distance += sqrt( pow( randomWalkCartesian[segment].x-intersection.x,2 ) + pow( randomWalkCartesian[segment].y-intersection.y,2 ) );
			printf("%.2f\t%.2f\tsolved\n", intersection.x, intersection.y);
			break;
		} else { // still in the forest
			printf("%.2f\t%.2f\tsolved\n",randomWalkCartesian[segment+1].x,randomWalkCartesian[segment+1].y);
			distance += randomWalkPolar[segment].r;
		}
	}
	translateCartesian(-1*points[worstePoint].x,-1*points[worstePoint].y);
	rotateCartesian(-1 * worsteAngle * 2 * M_PI / numberTestAngles);
	return distance;
}

/////////////////////////////////////////////////////////////////////////////////////
// genetic algorithm functions
/////////////////////////////////////////////////////////////////////////////////////

void makePopulation() {
	for (int i = 0; i < populationSize; i++) {
		makeWalk();
        for (int j = 0; j < degreesOfFreedom; j++)
            population[i][j] = randomWalkPolar[j];
    }
}

void immigration() {
	for (int immigrint = children + survivors; immigrint < populationSize; immigrint++) {
		makeWalk();
        for (int j = 0; j < degreesOfFreedom; j++)
            newPopulation[immigrint][j] = randomWalkPolar[j];
    }
}

void takeFromPopulation(int i) {
	for (int j = 0; j < degreesOfFreedom; j++ )
		randomWalkPolar[j] = population[i][j];
	polarToCartesian();
}

void populationFitness() {
	for (int i = 0; i < populationSize; i++) {
		takeFromPopulation(i);
		fitness[i] = worsteCaseScenario(&worsteAngle,&worstePoint);
	}
}

void findMinFitness() {
	minFitnessIndex = 0;
	for (int i = 0; i < populationSize; i++)
		if (fitness[i] < fitness[minFitnessIndex])
			minFitnessIndex = i;
}

int parentSelection() {
    // Calculate total inverse fitness
    double totalInverseFitness = 0.0;
    for (int i = 0; i < populationSize; i++)
        totalInverseFitness += 1.0 / ( fitness[i] );

    // Generate a random number between 0 and totalInverseFitness
    double randNum = ((double)rand() / RAND_MAX) * totalInverseFitness;

    // Find the index of the selected parent based on the cumulative inverse fitness
    double cumulativeInverseFitness = 0.0;
    for (int i = 0; i < populationSize; i++) {
        cumulativeInverseFitness += 1.0 / ( fitness[i] );
        if (cumulativeInverseFitness >= randNum) {
            return i; // Return the index of the selected parent
        }
    }

    // This line should not be reached, but in case of an issue, return the last index
    return populationSize - 1;
}

// insert a new x,y point in the list of cartesian coordinates. Once the new point
void insertCartesianAfter(int n, float x, float y) {
	for (int i = degreesOfFreedom; i > n + 1 ; i-- )
		randomWalkCartesian[i] = randomWalkCartesian[i-1];
	
	// insert the new point
	randomWalkCartesian[n + 1].x = x;
	randomWalkCartesian[n + 1].y = y;
	degreesOfFreedom++;
	
	// make sure the polar format matches and remove the last point
	cartesianToPolar();
	degreesOfFreedom--;
	randomWalkPolar[degreesOfFreedom-1].r = 300;
	polarToCartesian();
	}

void crossover(int parent1, int parent2, int child) {
    
	float parameterSetOne[200];
	float parameterSetTwo[200];
	float parameterSetNew[200];
	
	for (int i = 0; i < degreesOfFreedom * 2; i += 2) {
		parameterSetOne[i+0] = population[parent1][i/2].r;
		parameterSetOne[i+1] = population[parent1][i/2].theta;
		parameterSetTwo[i+0] = population[parent2][i/2].r;
		parameterSetTwo[i+1] = population[parent2][i/2].theta;
	}
	
	int crossoverPoint = rand() % ( degreesOfFreedom * 2 );
    for (int i = 0; i < crossoverPoint; i++)
        parameterSetNew[i] = parameterSetOne[i];
    for (int i = crossoverPoint; i < degreesOfFreedom * 2; i++)
        parameterSetNew[i] = parameterSetTwo[i];
  
	for (int i = 0; i < degreesOfFreedom; i++) {
		newPopulation[child][i].r = parameterSetNew[(i*2)+0];
		newPopulation[child][i].theta = parameterSetNew[(i*2)+1];
	}
}
	
void keepSurvivors() {
	for (int survivor = 0; survivor < survivors; survivor++) {
		findMinFitness();
		for (int i = 0; i < degreesOfFreedom; i++)
			newPopulation[children + survivor][i] = population[minFitnessIndex][i];
		fitness[minFitnessIndex] = 1000;
    }
}

void replacePopulation() {
	for (int i = 0; i < populationSize; i++)
        for (int j = 0; j < degreesOfFreedom; j++)
            population[i][j] = newPopulation[i][j];
}
	
void insertRandomPoint(int escapeSegment) {
    Cartesian p1 = randomWalkCartesian[escapeSegment];
    Cartesian p2 = randomWalkCartesian[escapeSegment + 1];
    Cartesian result;

    // Generate a random t value between 0 and 1
    float t = (float)rand() / RAND_MAX;

    // Interpolate between the two consecutive points using t
    result.x = p1.x + t * (p2.x - p1.x);
    result.y = p1.y + t * (p2.y - p1.y);
	result.x = result.x + (float)rand() / RAND_MAX;
	result.y = result.y + (float)rand() / RAND_MAX;
	insertCartesianAfter(escapeSegment,result.x,result.y);
}
	
void geneInsertion(int child) {
	for (int j = 0; j < degreesOfFreedom; j++ )
		randomWalkPolar[j] = newPopulation[child][j];
	polarToCartesian();
	worsteCaseScenario(&worsteAngle,&worstePoint);
	
	// get the escape segment and escape length for the path
	Cartesian intersection;
	int escapeSegment = 0;
	float escDistance = 0;
	rotateCartesian(worsteAngle * 2 * M_PI / numberTestAngles);
	translateCartesian(points[worstePoint].x,points[worstePoint].y);
	for (escapeSegment = 0; escapeSegment < degreesOfFreedom; escapeSegment++) {
		if ( detectEscape(randomWalkCartesian[escapeSegment],randomWalkCartesian[escapeSegment+1],&intersection) ) {
			escDistance += sqrt( pow( randomWalkCartesian[escapeSegment].x-intersection.x,2 ) + pow( randomWalkCartesian[escapeSegment].y-intersection.y,2 ) );
			break;
		} else { // still in the forest
			escDistance += randomWalkPolar[escapeSegment].r;
		}
	}
	polarToCartesian();

	// compute an insert distance and find which segment it's on
	int insertSegment = 0;
	float insDistance = 0;
	insDistance = escDistance * (float)rand() / RAND_MAX;
	for (insertSegment = 0; insertSegment < degreesOfFreedom; insertSegment++) {
		insDistance = insDistance - sqrt( pow( randomWalkCartesian[insertSegment].x-randomWalkCartesian[insertSegment+1].x,2 ) 
	                                    + pow( randomWalkCartesian[insertSegment].y-randomWalkCartesian[insertSegment+1].y,2 ) );
		if ( insDistance < 0 )
			break;
	}
	
	// insert a random point on the insertion segment
	insertRandomPoint(insertSegment);
	for (int j = 0; j < degreesOfFreedom; j++ )
		newPopulation[child][j] = randomWalkPolar[j];
}

void geneMutations(int child) {
	for (int j = 0; j < degreesOfFreedom; j++ )
		randomWalkPolar[j] = newPopulation[child][j];
	polarToCartesian();
	for (int j = 1; j < degreesOfFreedom; j++ ) {
		if ( rand() % 2 ) {
			randomWalkCartesian[j].x += (float)rand() / RAND_MAX;
			randomWalkCartesian[j].y += (float)rand() / RAND_MAX;
		}
	}
	cartesianToPolar();
	for (int j = 0; j < degreesOfFreedom; j++ )
		newPopulation[child][j] = randomWalkPolar[j];
	}

/////////////////////////////////////////////////////////////////////////////////////
// main programming section
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

	// Seed the random number generator with the current time
	int seed = time(NULL);
    srand(seed);
	
	// generate a default forest (user may overwrite with --forest or --hull options)
	makeForest();
	
	////////////////////////////////////////////////////////////////////////////
	// get the parameters of the run from the command line
	////////////////////////////////////////////////////////////////////////////
	int c;
    int option_index    = 0;
    struct option long_options[] = {
		{"segments",   required_argument, 0, 'j' },
		{"population", required_argument, 0, 'p' },
		{"angles",     required_argument, 0, 'a' },
		{"generations",required_argument, 0, 't' },
		{"children",   required_argument, 0, 'c' },
		{"gridsize",   required_argument, 0, 'g' },
		{"forest",     required_argument, 0, 'f' },
		{"hull",       required_argument, 0, 'h' },
		{"survivors",  required_argument, 0, 's' },
		{0,            0,                 0, 0   }
	};

	while (1) {
		   c = getopt_long(argc, argv,"",long_options, &option_index);
		   if (c == -1) 
			   break;
		   switch (c) {
			    case 0:
					printf("option %s", long_options[option_index].name);
					if (optarg)
						printf(" with arg %s", optarg);
					printf("\n");
					break;
				case 'a':
					numberTestAngles = atoi(optarg);
					break;
				case 'j':
					degreesOfFreedom = atoi(optarg);
					break;
				case 'p':
					populationSize = atoi(optarg);
					break;
				case 't':
					generations = atoi(optarg);
					break;
				case 'c':
					children = atoi(optarg);
					break;
				case 'g':
					forestGridSize = atof(optarg);
					break;	
				case 'f':
					readPointsFromFile(optarg);
					break;
				case 's':
					survivors = atoi(optarg);
					break;
				case 'h':
					totalForestSides = atoi(optarg);
					makeForest();
					break;
		   		case '?':
					break;
				default:
					printf("?? getopt returned character code 0%o ??\n", c);
					break;
               }
           }
		   
	////////////////////////////////////////////////////////////////////////////
	// quick check for possible input errors
	//////////////////////////////////////////////////////////////////////////// 
	int errors = 0;
	if ( totalForestSides > 64 )
		printf("error %d, forest size must be less than 64\n",++errors);
	if ( populationSize > 512 )
		printf("error %d, population size must be less than 512\n",++errors);
	if ( degreesOfFreedom > 64 )
		printf("error %d, escape strategy must have less than 64 sides\n",++errors);
	if ( children > populationSize )
		printf("error %d, offspring must not exceed the population\n",++errors);
	if ( ( children + survivors ) > populationSize )
		printf("error %d, children and offspring must not exceed the population\n",++errors);
	if ( errors )
		return 0;
		   
	////////////////////////////////////////////////////////////////////////////
	// run the genetic algorithm
	//////////////////////////////////////////////////////////////////////////// 	 
	
	// make test points in the forest
	gridForest();
	printForest();
	printGrid();
	
	// use genetic algorithm to find best escape path
	makePopulation();
	
	for (int generation = 0; generation < generations; generation++) {
        populationFitness();
		findMinFitness();
		fprintf(stderr, "generation = %d, bestIndex = %d, bestFitness = %f\n",generation,minFitnessIndex,fitness[minFitnessIndex]);
		printf("%.2f\t%.2f\tgeneration\n", (float)generation, (float)fitness[minFitnessIndex]);
		
		for (int child = 0; child < children; child = child + 2) {
			int parent1 = parentSelection(fitness);
			int parent2 = parentSelection(fitness);
			crossover(parent1, parent2, child + 0);
			crossover(parent2, parent1, child + 1);
			geneInsertion(child + 0);
			geneMutations(child + 1);
		}
		
	// housekeeping for the next generation
	keepSurvivors();
	immigration();
	replacePopulation();
	}
	
	populationFitness();
	findMinFitness();
	takeFromPopulation(minFitnessIndex);
	printCartesian();
	printRotations();
	float worsteCase = worsteCaseScenario(&worsteAngle,&worstePoint);
	float check = printSolution(worsteAngle,worstePoint);
	printf("%.2f\t%.2f\tsystem\n", (float)worsteCase, (float)degreesOfFreedom);
	fprintf(stderr,"final distance %.3f, check distance %.3f\n",worsteCase,check);
									   
return 0;
}




