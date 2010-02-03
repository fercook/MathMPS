/*
 * This file automatically produced by /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions/mprep from:
 *	mathematicatemplate.tm
 * mprep Revision 13 Copyright (c) Wolfram Research, Inc. 1990-2008
 */

#define MPREP_REVISION 13

#include "mathlink.h"

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

/********************************* end header *********************************/


# line 1 "mathematicatemplate.tm"
void findgroundMPSsite (int ,int ,int ,int ,int );

# line 32 "mathematicatemplatetm.c"


# line 13 "mathematicatemplate.tm"
int isLinkActive (void);

# line 38 "mathematicatemplatetm.c"


# line 25 "mathematicatemplate.tm"
void matrixTimesVector(int ,int ,int ,int ,int );
# line 43 "mathematicatemplatetm.c"


# line 34 "mathematicatemplate.tm"
void testCommunication(int ,int ,int ,int ,int );

# line 49 "mathematicatemplatetm.c"


void findgroundMPSsite P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5));

#if MLPROTOTYPES
static int _tr0( MLINK mlp)
#else
static int _tr0(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;

	findgroundMPSsite(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
 L4: L3: L2: L1: 
L0:	return res;
} /* _tr0 */


int isLinkActive P(( void));

#if MLPROTOTYPES
static int _tr1( MLINK mlp)
#else
static int _tr1(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp0;
	if ( ! MLNewPacket(mlp) ) goto L0;

	_tp0 = isLinkActive();

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutInteger( mlp, _tp0);

L0:	return res;
} /* _tr1 */


void matrixTimesVector P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5));

#if MLPROTOTYPES
static int _tr2( MLINK mlp)
#else
static int _tr2(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;

	matrixTimesVector(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
 L4: L3: L2: L1: 
L0:	return res;
} /* _tr2 */


void testCommunication P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5));

#if MLPROTOTYPES
static int _tr3( MLINK mlp)
#else
static int _tr3(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;

	testCommunication(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
 L4: L3: L2: L1: 
L0:	return res;
} /* _tr3 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((MLINK));
	char  *f_name;
	} _tramps[4] = {
		{ 5, 2, _tr0, "findgroundMPSsite" },
		{ 0, 0, _tr1, "isLinkActive" },
		{ 5, 2, _tr2, "matrixTimesVector" },
		{ 5, 2, _tr3, "testCommunication" }
		};

static char* evalstrs[] = {
	"FindGroundMPSSite::usage = \"FindGroundMPSSite[MPSsite,DLeft,DRig",
	"ht,hLeft,hRight,vLeft,vRight,Hamiltonian] \"",
	(char*)0,
	"IsLinkActive::usage = \"IsLinkActive[param] returns 1 if the link",
	" is active and working, anything else if it is not\"",
	(char*)0,
	(char*)0
};
#define CARDOF_EVALSTRS 2

static int _definepattern P(( MLINK, char*, char*, int));

static int _doevalstr P(( MLINK, int));

int  _MLDoCallPacket P(( MLINK, struct func[], int));


#if MLPROTOTYPES
int MLInstall( MLINK mlp)
#else
int MLInstall(mlp) MLINK mlp;
#endif
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = _definepattern(mlp, "FindGroundMPSSite[cA_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]", "{Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]], Normal[Re[cA]],Normal[Im[cA]], Normal[Re[cDLeft]], Normal[Im[cDLeft]], Normal[Re[cDRight]], Normal[Im[cDRight]], Normal[Re[chLeft]], Normal[Im[chLeft]], Normal[Re[chRight]], Normal[Im[chRight]], Normal[Re[cvLeft]], Normal[Im[cvLeft]], Normal[Re[cvRight]], Normal[Im[cvRight]], Normal[cHam]}", 0);
	if (_res) _res = _doevalstr( mlp, 0);
	if (_res) _res = _definepattern(mlp, "IsLinkActive[]", "{}", 1);
	if (_res) _res = _doevalstr( mlp, 1);
	if (_res) _res = _definepattern(mlp, "MPSMatrixVector[cA_List,cOut_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]", "{Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]], Normal[Re[cA]],Normal[Im[cA]], Normal[Re[cOut]],Normal[Im[cOut]], Normal[Re[cDLeft]], Normal[Im[cDLeft]], Normal[Re[cDRight]], Normal[Im[cDRight]], Normal[Re[chLeft]], Normal[Im[chLeft]], Normal[Re[chRight]], Normal[Im[chRight]], Normal[Re[cvLeft]], Normal[Im[cvLeft]], Normal[Re[cvRight]], Normal[Im[cvRight]], Normal[cHam]}", 2);
	if (_res) _res = _definepattern(mlp, "TestCommunication[cA_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]", "{Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]],Normal[Re[cA]],Normal[Im[cA]]}", 3);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


#if MLPROTOTYPES
int MLDoCallPacket( MLINK mlp)
#else
int MLDoCallPacket( mlp) MLINK mlp;
#endif
{
	return _MLDoCallPacket( mlp, _tramps, 4);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif


#if CARDOF_EVALSTRS
#if MLPROTOTYPES
static int  _doevalstr( MLINK mlp, int n)
#else
static int  _doevalstr( mlp, n)
	 MLINK mlp; int n;
#endif
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, charsleft, bytesleft);
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


#if MLPROTOTYPES
static int  _definepattern( MLINK mlp, char *patt, char *args, int func_n)
#else
static int  _definepattern( mlp, patt, args, func_n)
	MLINK  mlp;
	char  *patt, *args;
	int    func_n;
#endif
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


#if MLPROTOTYPES
int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
#else
int _MLDoCallPacket( mlp, functable, nfuncs)
	MLINK mlp;
	struct func functable[];
	int nfuncs;
#endif
{
	long len;
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! MLCheckFunction(mlp, "List", &len)
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* _MLDoCallPacket */


#if MLPROTOTYPES
mlapi_packet MLAnswer( MLINK mlp)
#else
mlapi_packet MLAnswer( mlp)
	MLINK mlp;
#endif
{
	mlapi_packet pkt = 0;

	while( !MLDone && !MLError(mlp) && (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
	MLAbort = 0;
	return pkt;
} /* MLAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

#if MLPROTOTYPES
static int refuse_to_be_a_frontend( MLINK mlp)
#else
static int refuse_to_be_a_frontend( mlp)
	MLINK mlp;
#endif
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluate( MLINK mlp, char *s)
#else
int MLEvaluate( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluate( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
} /* MLEvaluate */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluateString( MLINK mlp, char *s)
#else
int MLEvaluateString( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluateString( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


#if MLINTERFACE >= 3
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, int message, int n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	int message, n;
#endif
#else
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, unsigned long message, unsigned long n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	unsigned long message, n;
#endif
#endif /* MLINTERFACE >= 3 */
{
	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}

#if MLPROTOTYPES
#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
static int _MLMain( argv, argv_end, commandline)
#if MLINTERFACE >= 3
  char **argv, argv_end;
  char *commandline;
#else
  charpp_ct argv, argv_end;
  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
	if( stdenv == (MLEnvironment)0) goto R0;

#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv, MLDefaultHandler, 0);
#endif /* MLINTERFACE >= 3 */


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE >= 3
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#else
		: MLOpenArgv( stdenv, argv, argv_end, &err);
#endif
	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
	if( stdhandler) MLSetMessageHandler( mlp, stdhandler);

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMainString( commandline)  char *commandline;
#else
int MLMainString( commandline)  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
}

#if MLPROTOTYPES
int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
#else
int MLMainArgv( argv, argv_end)  char **argv, **argv_end;
#endif
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);

}

#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMain( int argc, char **argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMain( argc, argv) int argc; char **argv;
#else
int MLMain( argc, argv) int argc; charpp_ct argv;
#endif /* MLINTERFACE >= 3 */
#endif
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
 
