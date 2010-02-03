/* To launch this program from within Mathematica use:
 *   In[1]:= link = Install["arpackformps"]
 *
 */

#include "mathlink.h"

extern int isLinkActive( int param);

extern void findgroundMPSsite (int spin,int sizeL,int sizeR,int rangeL,int rangeR);

extern void matrixTimesVector(int spin,int sizeL,int sizeR,int rangeL,int rangeR);

extern void testCommunication (int spin,int sizeL,int sizeR,int rangeL,int rangeR);

extern "C"
{
  // Fortran subroutine:
  // ComputeGroundState(Ar,Ai,Leftr,Lefti,Rightr,Righti,DLeftr,Dlefti,DRightr,DRighti, &
  //         & hLeftr,hLefti,hRightr,hRighti,vLeftr,vLefti,vRightr,vRighti,Ham,spin,sizeL,sizeR,rLeft,rRight,MyInfo)
  void arpackdriver_mp_computegroundstate_(double * Ar, double * Ai,
					    double * DLeftr, double * DLefti, double * DRightr, double * DRighti,
					    double * hLeftr, double * hLefti, double * hRightr, double * hRighti,
					    double * vLeftr, double * vLefti, double * vRightr, double * vRighti,
					    double *Ham, int* spin,int* sizeL,int* sizeR,int* rLeft,int* rRight,int* MyInfo, double * energy);
}


extern "C"
{
  // Fortran subroutine:
  // ComputeGroundState(Ar,Ai,Leftr,Lefti,Rightr,Righti,DLeftr,Dlefti,DRightr,DRighti, &
  //         & hLeftr,hLefti,hRightr,hRighti,vLeftr,vLefti,vRightr,vRighti,Ham,spin,sizeL,sizeR,rLeft,rRight,MyInfo)
    void arpackdriver_mp_matrixtimesvectorfromc_(double * inputr, double * inputi ,double * outputr, double * outputi,
                                            double * DLeftr, double * DLefti, double * DRightr, double * DRighti,
                                            double * hLeftr, double * hLefti, double * hRightr, double * hRighti,
                                            double * vLeftr, double * vLefti, double * vRightr, double * vRighti,
                                            double *Ham, int* spin,int * vectorsize, int* sizeL,int* sizeR,int* rLeft,int* rRight);
}
// This is the wrapper for the (F90) main computational routine
// It gets data from Mathematica, sends it to the Fortran routine, and sends 
// results back to Mathematica.
void findgroundMPSsite(int spin,int sizeL,int sizeR,int rangeL,int rangeR)
{
  //Remaining arguments to get manually: Re[A], Im[A], Re[Left], Im[Left], Re[Right], Im[Right], Re[DLeft], Im[Dleft], 
  //                                     Re[DRight], Im[DRight], Re[hLeft], Im[hLeft], Re[hRight], Im[hRight], 
  //                                     Re[vLeft], Im[vLeft], Re[vRight], Im[vRight], Ham
  //


  // I will need to declare all the dimensional variables of all the input so I can free their memory later
  
  int *dimsAr,*dimsAi;
  int *dimsDLeftr,*dimsDLefti,*dimsDRightr,*dimsDRighti;
  int *dimshLeftr,*dimshLefti,*dimshRightr,*dimshRighti;
  int *dimsvLeftr,*dimsvLefti,*dimsvRightr,*dimsvRighti,*dimsHam;

  char **headsAr,**headsAi;
  char **headsDLeftr,**headsDLefti,**headsDRightr,**headsDRighti;
  char **headshLeftr,**headshLefti,**headshRightr,**headshRighti;
  char **headsvLeftr,**headsvLefti,**headsvRightr,**headsvRighti,**headsHam;

  int depthAr,depthAi;
  int depthDLeftr,depthDLefti,depthDRightr,depthDRighti;
  int depthhLeftr,depthhLefti,depthhRightr,depthhRighti;
  int depthvLeftr,depthvLefti,depthvRightr,depthvRighti,depthHam;

  double *Ar,*Ai;
  double *DLeftr,*DLefti,*DRightr,*DRighti;
  double *hLeftr,*hLefti,*hRightr,*hRighti;
  double *vLeftr,*vLefti,*vRightr,*vRighti,*Ham;

  double energy;
  int n,info;

  // Start getting data
  // Re[A], Im[A]
  MLGetReal64Array(stdlink, &Ar, &dimsAr, &headsAr, &depthAr);
  MLGetReal64Array(stdlink, &Ai, &dimsAi, &headsAi, &depthAi);

  // Re[DLeft], Im[Dleft], 
  MLGetReal64Array(stdlink, &DLeftr, &dimsDLeftr, &headsDLeftr, &depthDLeftr);
  MLGetReal64Array(stdlink, &DLefti, &dimsDLefti, &headsDLefti, &depthDLefti);

  // Re[DRight], Im[DRight]
  MLGetReal64Array(stdlink, &DRightr, &dimsDRightr, &headsDRightr, &depthDRightr);
  MLGetReal64Array(stdlink, &DRighti, &dimsDRighti, &headsDRighti, &depthDRighti);

  // Re[hLeft], Im[hleft], 
  MLGetReal64Array(stdlink, &hLeftr, &dimshLeftr, &headshLeftr, &depthhLeftr);
  MLGetReal64Array(stdlink, &hLefti, &dimshLefti, &headshLefti, &depthhLefti);

  // Re[hRight], Im[hRight]
  MLGetReal64Array(stdlink, &hRightr, &dimshRightr, &headshRightr, &depthhRightr);
  MLGetReal64Array(stdlink, &hRighti, &dimshRighti, &headshRighti, &depthhRighti);
 
  // Re[vLeft], Im[vleft], 
  MLGetReal64Array(stdlink, &vLeftr, &dimsvLeftr, &headsvLeftr, &depthvLeftr);
  MLGetReal64Array(stdlink, &vLefti, &dimsvLefti, &headsvLefti, &depthvLefti);

  // Re[vRight], Im[vRight]
  MLGetReal64Array(stdlink, &vRightr, &dimsvRightr, &headsvRightr, &depthvRightr);
  MLGetReal64Array(stdlink, &vRighti, &dimsvRighti, &headsvRighti, &depthvRighti);

  // Ham
  MLGetReal64Array(stdlink, &Ham, &dimsHam, &headsHam, &depthHam);

  // With all the data ready, now call the Fortran subroutine to use ARPACK  

    arpackdriver_mp_computegroundstate_(Ar,Ai,DLeftr,DLefti,DRightr,DRighti,
				       hLeftr, hLefti, hRightr, hRighti,vLeftr, vLefti, vRightr, vRighti,
				       Ham, &spin, &sizeL, &sizeR, &rangeL, &rangeR, &info, &energy);
  
  // ARPACK returned the result in the Ar and Ai variables, return this converted into a complex tensor
    MLPutFunction(stdlink,"List",3);
	MLPutReal64(stdlink,energy);
    MLPutFunction(stdlink,"Plus",2);
    MLPutReal64Array(stdlink,Ar,dimsAr,headsAr,depthAr);
    MLPutFunction(stdlink,"Times",2);
    MLPutSymbol(stdlink,"I");
    MLPutReal64Array(stdlink,Ai,dimsAi,headsAi,depthAi);
	MLPutInteger(stdlink,info);

  // Finally, release all that memory
  //A
  MLReleaseReal64Array(stdlink,Ar,dimsAr,headsAr,depthAr);
  MLReleaseReal64Array(stdlink,Ai,dimsAi,headsAi,depthAi);
  //DLeft, DRight
  MLReleaseReal64Array(stdlink,DLeftr,dimsDLeftr,headsDLeftr,depthDLeftr);
  MLReleaseReal64Array(stdlink,DLefti,dimsDLefti,headsDLefti,depthDLefti);
  MLReleaseReal64Array(stdlink,DRightr,dimsDRightr,headsDRightr,depthDRightr);
  MLReleaseReal64Array(stdlink,DRighti,dimsDRighti,headsDRighti,depthDRighti);
  //hLeft, hRight
  MLReleaseReal64Array(stdlink,hLeftr,dimshLeftr,headshLeftr,depthhLeftr);
  MLReleaseReal64Array(stdlink,hLefti,dimshLefti,headshLefti,depthhLefti);
  MLReleaseReal64Array(stdlink,hRightr,dimshRightr,headshRightr,depthhRightr);
  MLReleaseReal64Array(stdlink,hRighti,dimshRighti,headshRighti,depthhRighti);
  //vLeft, vRight
  MLReleaseReal64Array(stdlink,vLeftr,dimsvLeftr,headsvLeftr,depthvLeftr);
  MLReleaseReal64Array(stdlink,vLefti,dimsvLefti,headsvLefti,depthvLefti);
  MLReleaseReal64Array(stdlink,vRightr,dimsvRightr,headsvRightr,depthvRightr);
  MLReleaseReal64Array(stdlink,vRighti,dimsvRighti,headsvRighti,depthvRighti);
  //Ham
  MLReleaseReal64Array(stdlink,Ham,dimsHam,headsHam,depthHam);
  
}


//This little routine is just used to check if the link is active and working.
int isLinkActive(void)
{
  return 1;
}

void matrixTimesVector(int spin,int sizeL,int sizeR,int rangeL,int rangeR)
{
  //Remaining arguments to get manually: Re[A], Im[A], Re[Left], Im[Left], Re[Right], Im[Right], Re[DLeft], Im[Dleft], 
  //                                     Re[DRight], Im[DRight], Re[hLeft], Im[hLeft], Re[hRight], Im[hRight], 
  //                                     Re[vLeft], Im[vLeft], Re[vRight], Im[vRight], Ham
    //


  // I will need to declare all the dimensional variables of all the input so I can free their memory later
  
    int *dimsAr,*dimsAi,*dimsOutr,*dimsOuti;
    int *dimsDLeftr,*dimsDLefti,*dimsDRightr,*dimsDRighti;
    int *dimshLeftr,*dimshLefti,*dimshRightr,*dimshRighti;
    int *dimsvLeftr,*dimsvLefti,*dimsvRightr,*dimsvRighti,*dimsHam;

    char **headsAr,**headsAi,**headsOutr,**headsOuti;
    char **headsDLeftr,**headsDLefti,**headsDRightr,**headsDRighti;
    char **headshLeftr,**headshLefti,**headshRightr,**headshRighti;
    char **headsvLeftr,**headsvLefti,**headsvRightr,**headsvRighti,**headsHam;

    int depthAr,depthAi,depthOutr,depthOuti;
    int depthDLeftr,depthDLefti,depthDRightr,depthDRighti;
    int depthhLeftr,depthhLefti,depthhRightr,depthhRighti;
    int depthvLeftr,depthvLefti,depthvRightr,depthvRighti,depthHam;

    double *Ar,*Ai,*Outr,*Outi;
    double *DLeftr,*DLefti,*DRightr,*DRighti;
    double *hLeftr,*hLefti,*hRightr,*hRighti;
    double *vLeftr,*vLefti,*vRightr,*vRighti,*Ham;

    int n,info,vectorsize;

  // Start getting data
  // Re[A], Im[A]
    MLGetReal64Array(stdlink, &Ar, &dimsAr, &headsAr, &depthAr);
    MLGetReal64Array(stdlink, &Ai, &dimsAi, &headsAi, &depthAi);

      // Start getting data
  // Re[Out], Im[Out]
    MLGetReal64Array(stdlink, &Outr, &dimsOutr, &headsOutr, &depthOutr);
    MLGetReal64Array(stdlink, &Outi, &dimsOuti, &headsOuti, &depthOuti);

  // Re[DLeft], Im[Dleft], 
    MLGetReal64Array(stdlink, &DLeftr, &dimsDLeftr, &headsDLeftr, &depthDLeftr);
    MLGetReal64Array(stdlink, &DLefti, &dimsDLefti, &headsDLefti, &depthDLefti);

  // Re[DRight], Im[DRight]
    MLGetReal64Array(stdlink, &DRightr, &dimsDRightr, &headsDRightr, &depthDRightr);
    MLGetReal64Array(stdlink, &DRighti, &dimsDRighti, &headsDRighti, &depthDRighti);

  // Re[hLeft], Im[hleft], 
    MLGetReal64Array(stdlink, &hLeftr, &dimshLeftr, &headshLeftr, &depthhLeftr);
    MLGetReal64Array(stdlink, &hLefti, &dimshLefti, &headshLefti, &depthhLefti);

  // Re[hRight], Im[hRight]
    MLGetReal64Array(stdlink, &hRightr, &dimshRightr, &headshRightr, &depthhRightr);
    MLGetReal64Array(stdlink, &hRighti, &dimshRighti, &headshRighti, &depthhRighti);
 
  // Re[vLeft], Im[vleft], 
    MLGetReal64Array(stdlink, &vLeftr, &dimsvLeftr, &headsvLeftr, &depthvLeftr);
    MLGetReal64Array(stdlink, &vLefti, &dimsvLefti, &headsvLefti, &depthvLefti);

  // Re[vRight], Im[vRight]
    MLGetReal64Array(stdlink, &vRightr, &dimsvRightr, &headsvRightr, &depthvRightr);
    MLGetReal64Array(stdlink, &vRighti, &dimsvRighti, &headsvRighti, &depthvRighti);

  // Ham
    MLGetReal64Array(stdlink, &Ham, &dimsHam, &headsHam, &depthHam);

    vectorsize=spin*sizeL*sizeR;
    arpackdriver_mp_matrixtimesvectorfromc_(Ar, Ai, Outr, Outi,DLeftr, DLefti, DRightr, DRighti,
	 hLeftr, hLefti, hRightr, hRighti, vLeftr, vLefti, vRightr, vRighti,
         Ham, &spin, &vectorsize, &sizeL, &sizeR, &rangeL, &rangeR);
  
  // ARPACK returned the result in the Ar and Ai variables, return this converted into a complex tensor
    MLPutFunction(stdlink,"Plus",2);
    MLPutReal64Array(stdlink,Outr,dimsOutr,headsOutr,depthOutr);
    MLPutFunction(stdlink,"Times",2);
    MLPutSymbol(stdlink,"I");
    MLPutReal64Array(stdlink,Outi,dimsOuti,headsOuti,depthOuti);

  // Finally, release all that memory
  //A
    MLReleaseReal64Array(stdlink,Ar,dimsAr,headsAr,depthAr);
    MLReleaseReal64Array(stdlink,Ai,dimsAi,headsAi,depthAi);
  //A
    MLReleaseReal64Array(stdlink,Outr,dimsOutr,headsOutr,depthOutr);
    MLReleaseReal64Array(stdlink,Outi,dimsOuti,headsOuti,depthOuti);
  //DLeft, DRight
    MLReleaseReal64Array(stdlink,DLeftr,dimsDLeftr,headsDLeftr,depthDLeftr);
    MLReleaseReal64Array(stdlink,DLefti,dimsDLefti,headsDLefti,depthDLefti);
    MLReleaseReal64Array(stdlink,DRightr,dimsDRightr,headsDRightr,depthDRightr);
    MLReleaseReal64Array(stdlink,DRighti,dimsDRighti,headsDRighti,depthDRighti);
  //hLeft, hRight
    MLReleaseReal64Array(stdlink,hLeftr,dimshLeftr,headshLeftr,depthhLeftr);
    MLReleaseReal64Array(stdlink,hLefti,dimshLefti,headshLefti,depthhLefti);
    MLReleaseReal64Array(stdlink,hRightr,dimshRightr,headshRightr,depthhRightr);
    MLReleaseReal64Array(stdlink,hRighti,dimshRighti,headshRighti,depthhRighti);
  //vLeft, vRight
    MLReleaseReal64Array(stdlink,vLeftr,dimsvLeftr,headsvLeftr,depthvLeftr);
    MLReleaseReal64Array(stdlink,vLefti,dimsvLefti,headsvLefti,depthvLefti);
    MLReleaseReal64Array(stdlink,vRightr,dimsvRightr,headsvRightr,depthvRightr);
    MLReleaseReal64Array(stdlink,vRighti,dimsvRighti,headsvRighti,depthvRighti);
  //Ham
    MLReleaseReal64Array(stdlink,Ham,dimsHam,headsHam,depthHam);
  
}


//This is only for testing purposes
void testCommunication(int spin,int sizeL,int sizeR,int rangeL,int rangeR)
{
    
    int *dimsAr,*dimsAi;
    char **headsAr,**headsAi;
    int depthAr,depthAi;
    double *Ar,*Ai;
		
    int n,info;

  // Start getting data
  // Re[A], Im[A]
    MLGetReal64Array(stdlink, &Ar, &dimsAr, &headsAr, &depthAr);
    MLGetReal64Array(stdlink, &Ai, &dimsAi, &headsAi, &depthAi);
  
    
    MLPutFunction(stdlink,"List",5);
    MLPutInteger(stdlink,spin);
    MLPutInteger(stdlink,sizeL);
    MLPutInteger(stdlink,sizeR);
    MLPutInteger(stdlink,rangeL);
    MLPutInteger(stdlink,rangeR);
    
    /*
    MLPutFunction(stdlink,"Plus",2);
    MLPutReal64Array(stdlink,Ar,dimsAr,headsAr,depthAr);
    MLPutFunction(stdlink,"Times",2);
    MLPutSymbol(stdlink,"2");
    MLPutReal64Array(stdlink,Ai,dimsAi,headsAi,depthAi);

*/
  // Finally, release all that memory
  //A
    MLReleaseReal64Array(stdlink,Ar,dimsAr,headsAr,depthAr);
    MLReleaseReal64Array(stdlink,Ai,dimsAi,headsAi,depthAi);

}


#if MACINTOSH_MATHLINK

int main( int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
	argc = argc; /* suppress warning */
	return MLMain( 0, argv);
}

#elif WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

	hinstPrevious = hinstPrevious; /* suppress warning */

	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( (int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

#endif
