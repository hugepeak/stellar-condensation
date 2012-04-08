#include <ctype.h>
#include <Libnucnet.h>
#include <WnSparseSolve.h>

#define D_CUTOFF    0.e-15
#define S_DEBUG     "debug"
#define I_BUF       256

int initialize_( const char *, const char * );
int getzone_( const char * );
int getnumberofspecies_( int * );
int printnumberofreactions_( void );
int getabundances_( int *, int *, double * );
int getelementalabundances_( double * );
int cleanup_( void );
int getabunds( Libnucnet__Zone *, void * );
int getspecies( Libnucnet__Species *, void * );
int decay_( double * );
int
evolveTimeStep(
  WnSparseSolve__Exp *, Libnucnet__Zone *, WnMatrix *, double
);
int zone_decay( Libnucnet__Zone *, void * );
int
clean_species( Libnucnet__Species *, Libnucnet__Zone * );
char * trimwhitespace( char * );

Libnucnet *p_libnucnet;

int
initialize_(
  const char *s_file_name,
  const char *s_nuc_xpath
)
{

   char *s1, *s2;

   p_libnucnet = Libnucnet__new();

   s1 = (char *) malloc( I_BUF * sizeof( char ) );
   s2 = (char *) malloc( I_BUF * sizeof( char ) );

   strcpy( s1, s_file_name );
   strcpy( s2, s_nuc_xpath );

   Libnucnet__Net__updateFromXml(
     Libnucnet__getNet( p_libnucnet ),
     trimwhitespace( s1 ),
     trimwhitespace( s2 ),
     NULL
   );

   free( s1 );
   free( s2 );

   return 0;
}

int
getnumberofspecies_( int *i_species )
{
  *i_species =
    (int)
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_libnucnet ) )
    );

  return 0;

}

int
printnumberofreactions_( void )
{
   printf(
     "Number reactions = %lu\n",
     (unsigned long)Libnucnet__Reac__getNumberOfReactions(
       Libnucnet__Net__getReac( Libnucnet__getNet( p_libnucnet ) )
     )
   ); 

   return 0;
}

int
cleanup_( void )
{
   Libnucnet__free( p_libnucnet );
   return 0;
}

int
getzone_( const char *s_zone_xml )
{

  char *s1;

  s1 = (char *) malloc( I_BUF * sizeof( char ) );

  strcpy( s1, s_zone_xml );

  Libnucnet__assignZoneDataFromXml(
    p_libnucnet,
    trimwhitespace( s1 ),
    NULL
  );

  if( Libnucnet__getNumberOfZones( p_libnucnet ) == 0 )
  {
    fprintf( stderr, "No zones were assigned!\n" );
    exit( EXIT_FAILURE );
  }

  free( s1 );

  return 0;

}

int
getabundances_( int *i_z, int *i_a, double *d_y )
{

  typedef struct {
    int *aZ;
    int *aA;
    double *aY;
  } user_data;

  user_data *p_user_data;

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) exit(0);

  p_user_data->aZ = (int *) i_z;
  p_user_data->aA = (int *) i_a;
  p_user_data->aY = (double *) d_y;

  Libnucnet__iterateZones(
    p_libnucnet,
    (Libnucnet__Zone__iterateFunction) getabunds,
    p_user_data
  );

  i_z = p_user_data->aZ;
  i_a = p_user_data->aA;
  d_y = p_user_data->aY;

  free( p_user_data );

  return 0;

}

int
getabunds( Libnucnet__Zone *p_zone, void *p_data )
{

  typedef struct {
    int *aZ;
    int *aA;
    double *aY;
  } user_data;

  typedef struct {
    user_data *pUserData;
    Libnucnet__Zone *pZone;
  } new_user_data;

  new_user_data *p_new_user_data;

  p_new_user_data = ( new_user_data * ) malloc( sizeof( new_user_data ) );

  if( !p_new_user_data ) exit(0);

  p_new_user_data->pUserData = ( user_data * ) p_data;
  p_new_user_data->pZone = p_zone;

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_libnucnet )
    ),
    (Libnucnet__Species__iterateFunction) getspecies,
    p_new_user_data
  );

  free( p_new_user_data );

  return 1;

}

int
getelementalabundances_( double *d_y_elem )
{

  size_t i, i_species;
  int *a_z, *a_a;
  double *a_y;
  
  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_libnucnet ) )
    );

  a_z = (int *) malloc( sizeof( int ) * i_species );
  a_a = (int *) malloc( sizeof( int ) * i_species );
  a_y = (double *) malloc( sizeof( double ) * i_species );

  getabundances_( a_z, a_a, a_y );

  for( i = 0; i < i_species; i++ )
  {
    if( a_z[i] > 0 ) d_y_elem[a_z[i]-1] = 0.;
  }

  for( i = 0; i < i_species; i++ )
  {
    if( a_z[i] > 0 ) d_y_elem[a_z[i]-1] += a_y[i];
  }

  free( a_z );
  free( a_a );
  free( a_y );

  return 0;

}

int
getspecies( Libnucnet__Species *p_species, void *p_data )
{

  typedef struct {
    int *aZ;
    int *aA;
    double *aY;
  } user_data;

  typedef struct {
    user_data *pUserData;
    Libnucnet__Zone *pZone;
  } new_user_data;

  size_t i_index;
  new_user_data *p_new_user_data;

  p_new_user_data = ( new_user_data * ) p_data;

  i_index = Libnucnet__Species__getIndex( p_species );

  p_new_user_data->pUserData->aZ[i_index] =
    (int) Libnucnet__Species__getZ( p_species );
  p_new_user_data->pUserData->aA[i_index] =
    (int) Libnucnet__Species__getA( p_species );
  p_new_user_data->pUserData->aY[i_index] =
    Libnucnet__Zone__getSpeciesAbundance(
      p_new_user_data->pZone,
      p_species
    );

  return 1;

}

int
decay_( double *d_tend )
{

  typedef struct {
    Libnucnet *pLibnucnet;
    double dTend;
    WnSparseSolve__Exp *pMySolver;
  } user_data;

  user_data *p_user_data;

  /*==========================================================================
  // Get user data
  //========================================================================*/

  p_user_data = (user_data *) malloc( sizeof( user_data ) );

  p_user_data->dTend = *d_tend;

  /*==========================================================================
  // Get solver.
  //========================================================================*/

  p_user_data->pMySolver = WnSparseSolve__Exp__new();

  if( !strcmp( S_DEBUG, "debug" ) )
    WnSparseSolve__Exp__setDebug( p_user_data->pMySolver );

  WnSparseSolve__Exp__updateMaximumIterations(
    p_user_data->pMySolver,
    1000000
  );

  WnSparseSolve__Exp__updateTolerance(
    p_user_data->pMySolver,
    1.e-16
  );

  /*==========================================================================
  // Get decay network.
  //========================================================================*/

  p_user_data->pLibnucnet = p_libnucnet;

  /*==========================================================================
  // Iterate zones.
  //========================================================================*/

  Libnucnet__iterateZones(
    p_user_data->pLibnucnet,
    (Libnucnet__Zone__iterateFunction) zone_decay,
    p_user_data
  );

  /*==========================================================================
  // Clean up.
  //========================================================================*/

  WnSparseSolve__Exp__free( p_user_data->pMySolver );

  p_libnucnet = p_user_data->pLibnucnet;

  free( p_user_data );

  return 0;

}

int
zone_decay(
  Libnucnet__Zone *p_zone, void *p_data
)
{

  typedef struct {
    Libnucnet *pLibnucnet;
    double dTend;
    WnSparseSolve__Exp *pMySolver;
  } user_data;

  double d_xsum_0;
  WnMatrix *p_matrix;
  user_data *p_user_data;

  /*==========================================================================
  // Get data.
  //========================================================================*/

  p_user_data = ( user_data * ) p_data;

  /*==========================================================================
  // Get initial mass fraction.
  //========================================================================*/

  d_xsum_0 = Libnucnet__Zone__computeAMoment( p_zone, 1 );

  /*==========================================================================
  // Compute rates
  //========================================================================*/

  Libnucnet__Zone__computeRates( p_zone, 1., 1. );

  /*==========================================================================
  // Get the Jacobian Matrix.
  //========================================================================*/

  p_matrix = Libnucnet__Zone__computeJacobianMatrix( p_zone );

  WnMatrix__scaleMatrix( p_matrix, -1. );


  /*==========================================================================
  // Evolution.
  //========================================================================*/

  if(
    evolveTimeStep(
      p_user_data->pMySolver,
      p_zone,
      p_matrix,
      p_user_data->dTend
    )
  )
      return 0;

  if( !strcmp( S_DEBUG, "debug" ) )
    printf(
      "t = %e  xsum0 - xsum = %e\n",
      p_user_data->dTend,
      d_xsum_0 - Libnucnet__Zone__computeAMoment( p_zone, 1 )
    );
   
  WnMatrix__free( p_matrix );

  /*==========================================================================
  // Clean up abundances.
  //========================================================================*/

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_user_data->pLibnucnet )
    ),
    (Libnucnet__Species__iterateFunction) clean_species,
    p_zone
  );

  return 1;

}

int
clean_species( Libnucnet__Species *p_species, Libnucnet__Zone *p_zone )
{

  if( Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species ) < D_CUTOFF )
    Libnucnet__Zone__updateSpeciesAbundance( p_zone, p_species, 0. );

  return 1;

}

int
evolveTimeStep(
  WnSparseSolve__Exp *p_my_solver,
  Libnucnet__Zone *self,
  WnMatrix *p_matrix,
  double d_dt
) {

  size_t i_rows;
  gsl_vector *p_yold, *p_sol;

  /*============================================================================
  // Get number of rows.
  //==========================================================================*/

  i_rows = WnMatrix__getNumberOfRows( p_matrix );

  /*============================================================================
  // Store old abundances.
  //==========================================================================*/

  p_yold = Libnucnet__Zone__getAbundances( self );

  /*============================================================================
  // Solve matrix.
  //==========================================================================*/

  p_sol = WnSparseSolve__Exp__solve( p_my_solver, p_matrix, p_yold, d_dt );

  /*============================================================================
  // Check for solution.
  //==========================================================================*/

  if( !p_sol ) return 1;

  /*============================================================================
  // Update abundances.
  //==========================================================================*/

  Libnucnet__Zone__updateAbundances( self, p_sol );
  gsl_vector_sub( p_sol, p_yold );
  Libnucnet__Zone__updateAbundanceChanges( self, p_sol );

  /*============================================================================
  // Free allocated memory.
  //==========================================================================*/

  gsl_vector_free( p_yold );
  gsl_vector_free( p_sol );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return 0;

}

/**
 * Routine to trim the trailing white spaces from a string.
 */

char *trimwhitespace( char *str )
{

  char *end;

/* Trim leading space */

  while(isspace(*str)) str++;

  if(*str == 0)  /* All spaces? */
    return str;

/* Trim trailing space */

  end = str + strlen(str) - 1;
  while(end > str && isspace(*end))
  {
    end--;
  }

/* Write new terminator */

  *(end+1) = 0;

  return str;

}
