void findgroundMPSsite (int ,int ,int ,int ,int );

:Begin: 
:Function:          findgroundMPSsite
:Pattern:           FindGroundMPSSite[cA_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]
:Arguments:         {Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]],Normal[Re[cA]],Normal[Im[cA]], Normal[Re[cDLeft]], Normal[Im[cDLeft]], Normal[Re[cDRight]], Normal[Im[cDRight]], Normal[Re[chLeft]], Normal[Im[chLeft]], Normal[Re[chRight]], Normal[Im[chRight]], Normal[Re[cvLeft]], Normal[Im[cvLeft]], Normal[Re[cvRight]], Normal[Im[cvRight]], Normal[cHam]}
:ArgumentTypes:     {Integer,Integer,Integer,Integer,Integer,Manual}
:ReturnType:        Manual
:End:

:Evaluate: FindGroundMPSSite::usage = "FindGroundMPSSite[MPSsite,DLeft,DRight,hLeft,hRight,vLeft,vRight,Hamiltonian] "

int isLinkActive (void);

:Begin: 
:Function:          isLinkActive
:Pattern:           IsLinkActive[]
:Arguments:         {}
:ArgumentTypes:     {}
:ReturnType:        Integer
:End:

:Evaluate: IsLinkActive::usage = "IsLinkActive[param] returns 1 if the link is active and working, anything else if it is not"

void matrixTimesVector(int ,int ,int ,int ,int );
:Begin: 
:Function:          matrixTimesVector
:Pattern:           MPSMatrixVector[cA_List,cOut_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]
:Arguments:         {Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]], Normal[Re[cA]],Normal[Im[cA]], Normal[Re[cOut]],Normal[Im[cOut]], Normal[Re[cDLeft]], Normal[Im[cDLeft]], Normal[Re[cDRight]], Normal[Im[cDRight]], Normal[Re[chLeft]], Normal[Im[chLeft]], Normal[Re[chRight]], Normal[Im[chRight]], Normal[Re[cvLeft]], Normal[Im[cvLeft]], Normal[Re[cvRight]], Normal[Im[cvRight]], Normal[cHam]}
:ArgumentTypes:     {Integer,Integer,Integer,Integer,Integer,Manual}
:ReturnType:        Manual
:End:

void testCommunication(int ,int ,int ,int ,int );

:Begin: 
:Function:               testCommunication
:Pattern:                TestCommunication[cA_List,cDLeft_List,cDRight_List,chLeft_List,chRight_List,cvLeft_List,cvRight_List,cHam_List]
:Arguments:         {Length[cA], Length[chLeft[[1]]], Length[chRight[[1]]], Length[cvLeft[[1]]], Length[cvRight[[1]]],Normal[Re[cA]],Normal[Im[cA]]}
:ArgumentTypes: {Integer,Integer,Integer,Integer,Integer,Manual}
:ReturnType:        Manual
:End: