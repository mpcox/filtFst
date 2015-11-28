/*

filtFst.c

Takes ms-style input (real or simulated) and returns Fst values using equation #3
of Hudson, Slatkin and Maddison (1992).

Calculates Fst between two groups with the first n1 sequences from group one
and the last n2 sequences from group two. Missing data can be encoded as '?'.

Code originally written by Jeff Wall, University of California San Francisco (2006)
Modified by August Woerner, University of Arizona (2007)
Modified again by Murray Cox, Massey University (2015)

Compilation: gcc filtFst.c -o filtFst -O3 -lm -Wall

*/


#include <stdio.h>
#include <stdlib.h>

short calculateFst(char **list, char *segSiteString);

int n1, n2, n, min_sample_size, ss, a, b;
int SSMAX;

int main(int argc, char *argv[]){
    
    char **list;
    SSMAX = 5000000;
    
    int isMs = 0; // boolean == 1 if input is simulated ms data
    char *segSiteString; // used to parse the number of segregating sites out of stdin
    char temp;
    
    // attempts to distinguish between simulated and real data in ms format
    if (!scanf("%d ", &ss)){
        
        isMs = 1; // if not simulated ms data, try the real thing
        segSiteString = "segsites: %d\n";
        
        // look at the ms command line
        while(1){
            
            temp = getchar();
            if(temp == 'n'){
                
                fprintf(stderr, "Failed to find the -I flag in your ms command line!\n");
                return 3;
                
            }
            else if(temp == '-'){
                
                temp = getchar();
                if(temp == 'I')
                break;              
            }
                 
        }
        
        int numberOfPops;
        if(scanf("%d %d %d", &numberOfPops, &n1, &n2) != 3){
            
            fprintf(stderr, "Failed to properly parse the ms -I flag!\n");
            return 4;            
        }
        else if(numberOfPops != 2){
            
            fprintf(stderr, "filtFst only works with two populations at a time, not %d\n", numberOfPops);
            return 5; // ensure that only two populations are specified
        }
        
        n = n1 + n2;        
    }
    else{
        
        segSiteString = "%d ";
        
        if(argc != 3 && argc != 4){
            
            fprintf(stderr, "Usage:\nfiltFst n1 n2 (Sample Size Per Population)\n");
            return 1;
        }
        
        n1 = atoi(argv[1]);
        n2 = atoi(argv[2]);
        n = n1 + n2;
        
        if(argc == 4){
            
            min_sample_size = atoi(argv[3]);
            // min = 2 to avoid divide by zero
            if (min_sample_size < 2)
            min_sample_size = 2;
            
        }
        else
        min_sample_size = 2;
        
        if(n1 <= min_sample_size || n2 <= min_sample_size){
            
            fprintf(stderr,
            "Illegal population sizes specified!\nMust be greater than %d\n", min_sample_size);
            fprintf(stderr, "Usage:\nfiltFst n1 n2 (Sample Size Per Population)\n");
            
            return 1;
        }
    }
    
    // allocate some memory
    list = (char **) malloc ((unsigned) n * sizeof (char *));
    if(list == NULL){
        
        fprintf(stderr, "Failed to malloc in main!\n");
        return 2;
    }
    
    for(a=0; a<n; ++a)
    if( (list[a] = (char *) malloc ((unsigned) SSMAX * sizeof (char)) ) == NULL)
    fprintf(stderr, "Failed to malloc in main!\n");
    
    if (isMs){
        
        // this parses the ms file itself
        while((temp = getchar()) != EOF){
            
            while(temp != EOF && temp != '/') // looks for the '//' delimiter
            temp = getchar();
            
            if(temp == EOF) {
                
                fprintf(stderr, "Premature ending of stdin\n");
                return 2;
            }
            
            getchar(); // second '/' delimiter
            getchar(); // newline
            
            if(scanf(segSiteString, &ss) != 1){
                
                fprintf(stderr, "Failed to grab the number of segregating sites from the ms input!\n");
                return 3;
            }
            
            scanf("positions: ");
            calculateFst(list, segSiteString);
        }
    }
    else{
        
        // read the 'pseudoMs' (i.e., real data files a la Jeff Wall)
        do{
            calculateFst(list, segSiteString);
        }
        while(scanf(segSiteString, &ss) == 1);
    }
    
    return 0;
}


short calculateFst(char **list, char *segSiteString){
    
    // dynamically scale the array
    if(ss >= SSMAX){
        
        while(ss >= SSMAX)
        SSMAX *= 2; // double the array size til it's larger than the number of segregating sites
        
        for(a = 0; a < n; ++a){
            
            list[a] = realloc(list[a], SSMAX);
            if(list[a] == NULL) {
                fprintf(stderr, "Failed to realloc in calculateFST!\n");
                exit(10);
            }
        }
    }
    
    int i1, i2, j1, j2, piIn1, piIn2, piBet, N1, N2;
    double Hw, Hb;
    
    if(ss == 0) // zero seg sites
    printf("0.0\n");
    else{
        
        for(a=0; a<ss; ++a)
        scanf("%*f ");
        
        for(a=0; a<n; ++a) // load the matrix
        if(scanf ("%s\n", list[a]) != 1) {
            // error handling added by AW 4/03/07
            fprintf(stderr, "Oops. scanf failed! Perhaps you specified the wrong sample size?\n");
            exit(10);
        }
        
        // the calculation here is from the original code.
        Hw = Hb = 0.;
        for(a=0; a<ss; ++a){
            
            piIn1 = piIn2 = piBet = 0;
            i1=i2=j1=j2=0;
            for(b=0; b<n1; ++b)
            if(list[b][a] == '1')
            ++i1;
            else if(list[b][a] == '0')
            ++i2;
            
            for(b=n1; b<n; ++b)
            if(list[b][a] == '1')
            ++j1;
            else if(list[b][a] == '0')
            ++j2;
            
            piIn1 += i1 * i2;
            piIn2 += j1 * j2;
            piBet = i1 * j2 + j1 * i2;
            N1 = i1 + i2;
            N2 = j1 + j2;
            
            // avoid divide by zero! AW 4/03/07
            // plus, Fst gets funky when there's a lot of missing data
            if(N1 >= min_sample_size && N2 >= min_sample_size){
                
                Hw += piIn1 / (double) (N1*(N1-1.));
                Hw += piIn2 / (double) (N2*(N2-1.));
                Hb += piBet / (double) (N1*N2);
            }
            
            // printf("%d %d %d %d %dn", a, piIn1, piIn2, piBet, ss);
        }
        
        if(Hb > 0)
            printf("%lf\n", 1. - (Hw / Hb));
        else
            printf("nan\n");
    }
    
    return 1;
}
