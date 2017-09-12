//
//  Nbody.c
//  HW3
//
//  Created by Dan Wilson on 11/11/2016.
//  Copyright © 2016 Dan Wilson. All rights reserved.
//

//This program simulates the motion of N bodies in free space which interact with eachother via the gravitational force, and then calculates the relative energy loss of each body due to innacuracies in the simulation methods. The program easily scales for number of bodies, however does not scale with dimensions. ITEMS_PER_LINE and the sscanf() within readOrbits must be adjusted accordingly if the number of dimensions is changed from 3, otherwise the program will scale properly. Under the given constants and input files, the program will simulate the stable orbit of the Sun, Earth and Moon for one year, using a timestep of a quarter day. The readOrbits() function was adapted from Charles Williams' code available on the Exeter University electronic learning environment (19/11/16).

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#define MAX_FILE_LINE_SIZE 250
#define ITEMS_PER_LINE 8
#define MAX_BODIES 3
#define MAX_NAME_SIZE 32

//defines the size of time step for the system, where 86164.1 is the number of seconds in a day, and 31536000 is the number of seconds in a year.
#define TIME_STEP (86164.1 / 4)
#define TIME_MAX (31536000 / 1)

#define G_CONST 6.67408e-11

//defines the names for the input file, the command for executing gnuplot, the script written by the program for plotting orbits, and the output data file.
#define INPUT_FILE "sun_earth_moon.txt"
#define GNUPLOT_EXE "gnuplot"
#define GNUPLOT_SCRIPT "gnuplot.script"
#define GNUPLOT_DATA "OrbitData_output.txt"

//creates the enumeration 'coords' for ease of use, and also to define the number of dimensions of the system.
typedef enum Coords { x, y, z } coords;

//creates the structure that contains all the relevant information of a particular body that is required to study it's motion due to gravity.
typedef struct body {
    char name[MAX_NAME_SIZE];
    double mass;
    double r[sizeof(coords)];
    double v[sizeof(coords)];
    double a[sizeof(coords)];
    double Energy;
} Body;

//function to check if there is enough memory when a structure is allocated
static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Out of memory!\n");
        exit(1);
    }
    return p;
}

//function to read the initial positions, masses and velocities of the bodies to be simulated from INPUT_FILE.
static int readOrbits(Body *bodies) {
    char line[MAX_FILE_LINE_SIZE];
    char nameDummy[MAX_FILE_LINE_SIZE];
    FILE *finput = fopen( INPUT_FILE, "r" );
    if (!finput) {
        fprintf(stderr, "Error: Could not open file '%s'.\n",INPUT_FILE);
        exit(1);
    }

    int bodyN = 0; /* Number of bodies successfully found */
    while ( bodyN < MAX_BODIES && fgets(line, MAX_FILE_LINE_SIZE, finput) ) {
        if (line[0] != '#') {
            int nItemsScanned = sscanf(line,"%s %lg %lg %lg %lg %lg %lg %lg",
                                       nameDummy, &bodies[bodyN].mass,
                                       &bodies[bodyN].r[x], &bodies[bodyN].r[y], &bodies[bodyN].r[z],
                                       &bodies[bodyN].v[x], &bodies[bodyN].v[y], &bodies[bodyN].v[z]);
            if (nItemsScanned == ITEMS_PER_LINE) {
                strncpy(bodies[bodyN].name,nameDummy,MAX_NAME_SIZE);
                bodyN++;
            } else {
                fprintf(stderr, "Unknown format: %s\n",line);
            }
        }
    }

    fclose(finput);
    return bodyN;

}

//function to write the script GNUPLOT_SCRIPT that contains the parameters for plotting the graph. The axis ranges are not specified in order to scale with varied input files. The graph is a surface plot with lines and points.
static void writeScript(int bodyN, Body * bodies) {

    FILE *fscript = fopen (GNUPLOT_SCRIPT, "w");
    if (!fscript) {
        fprintf(stderr, "Error: Could not open file '%s'.\n",INPUT_FILE);
        exit(1);
    }
    fprintf(fscript,    "set size square\n"
                        "set xlabel \"Distance (m)\"\n"
                        "splot ");
    for(int i=0; i<bodyN;i++) {

        fprintf(fscript, "'%s' using %d:%d:%d title '%s' with linespoints,\\\n" ,GNUPLOT_DATA, 3*i+1, 3*i+2, 3*i+3, bodies[i].name);
    }

    fclose(fscript);


}


//VERLET: function to calculate the total acceleration of a body by summing up the accelerations due to the gravitational force between all bodies of the system. The accelerations are initialised to zero before each sum.
static void aNew(int bodyN, Body *bodies) {
    double R;
    for (int i=0; i < bodyN; i++) {
        for (int k=0; k < sizeof(coords); k++) {
        bodies[i].a[k] = 0.0;
        }
    }

    for (int i=0; i<bodyN ;i++) {
        for(int j=0; j<bodyN ;j++) {
            if(i!=j) {
                R = sqrt( pow((bodies[i].r[x] - bodies[j].r[x]),2) + pow((bodies[i].r[y] - bodies[j].r[y]),2) + pow((bodies[i].r[z] - bodies[j].r[z]),2));
                for(int k=0; k<sizeof(coords); k++) {
                    bodies[i].a[k] += -G_CONST * bodies[j].mass * ((bodies[i].r[k] - bodies[j].r[k] ) / pow(R,3));
                }
            }
        }
    }
}

//VERLET: function to calculate the new position of the body given its velocity and the set time step. The function also prints the positions to a file once evaluated.
static void rNew(int bodyN, Body *bodies, FILE *foutput) {

    for(int i=0; i < bodyN; i++){
        for(int k=0; k<sizeof(coords); k++) {
            bodies[i].r[k] += bodies[i].v[k]*TIME_STEP;
        }
        fprintf(foutput, "%.8g\t%.8g\t%.8g\t", bodies[i].r[x], bodies[i].r[y], bodies[i].r[z] );
    }
    fprintf(foutput, "\n");
}

//VERLET: function to calculate the velocity of the body after half a time step. Ths function is called twice for each time step, and is calculated as such to minimise the number of variables that need to be stored.
static void vhalfNew(int bodyN, Body *bodies) {

    for(int i=0; i<bodyN; i++) {
        for(int k=0; k<sizeof(coords); k++) {
            bodies[i].v[k] += 0.5*bodies[i].a[k]*TIME_STEP;
        }
    }
}

//VERLET: function to call the individual equations required for the verlet algorithm for time TIME_MAX in steps of TIME_STEP.
static void verlet(int bodyN, Body * bodies, FILE *foutput) {

    for (int k= 0; k < (TIME_MAX / TIME_STEP) ; k++) {

        vhalfNew(bodyN, bodies);
        rNew(bodyN, bodies, foutput);
        aNew(bodyN, bodies);
        vhalfNew(bodyN, bodies);

    }

}

//function to evaluate the total energy of each body
static void energyCheck(int bodyN, double *E, Body * bodies) {

    for(int i=0; i<bodyN; i++) {
        for(int j=0; j<bodyN; j++) {
            if(i!=j) {
                E[i] += - G_CONST * bodies[i].mass * bodies[j].mass / (2 * sqrt( pow((bodies[i].r[x] - bodies[j].r[x]),2) + pow((bodies[i].r[y] - bodies[j].r[y]),2) + pow((bodies[i].r[z] - bodies[j].r[z]),2)));
            }
        }
    }
}


int main(int argc, const char * argv[]) {


    Body *bodies = xmalloc(MAX_BODIES * sizeof(Body));
    int bodyN = readOrbits(bodies);


    double Eold[MAX_BODIES] = {0};
    energyCheck(bodyN, Eold, bodies);

    aNew(bodyN,bodies);

    FILE *foutput = fopen (GNUPLOT_DATA, "w");
    if (!foutput) {
        fprintf(stderr, "Error: Could not open file '%s'.\n",INPUT_FILE);
        exit(1);
    }


    verlet(bodyN, bodies, foutput);


    fclose(foutput);



    double Enew[MAX_BODIES] = {0};
    energyCheck(bodyN, Enew, bodies);

    printf("Body\tAbsolute Energy Loss\n");
    for( int i=0; i<bodyN; i++) {
        printf("%s\t%.2g\n", bodies[i].name, fabs((Eold[i] - Enew[i]) / Eold[i])  );
    }



    writeScript(bodyN, bodies);

    char command[PATH_MAX];
    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
    system( command );

    free(bodies);
}
