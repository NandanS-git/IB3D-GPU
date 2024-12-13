
      PROGRAM main
        USE global
        IMPLICIT NONE
        CALL readInput
        CALL readSurfaceMeshGmsh
        !CALL readSurfaceMeshGambit
        CALL allocateArrays
        CALL shiftSurfaceNodesInitial
        CALL computeSurfaceNorm
        CALL tagging
        CALL writeTagging
        CALL cellCount
        CALL computeNormDistance
        IF (iStart.eq.0) CALL initialConditions 
        IF (iStart.eq.1) CALL lastConditions       
        CALL coefficientMatrix
        CALL non_uni_coeff
        ita = 1
        ita1 = 1
        totime = 0.
        !CALL readComputeSumData
        !CALL readstressData
        !ita = ita + 1
        !ita1 = ita1 + 1
        totime = totime + deltat
        !CALL calculateUcMO
        !CALL calculateUc               
        CALL nsMomentum               
        CALL velocityBC        
        CALL solidCellBC        
        CALL velocityForcing1         
        CALL velocityBC        
        CALL poissonSolver        
        CALL pressureForcing1              
        !CALL calculateUcMO
        !CALL calculateUc      
        !CALL computeSumData
        CALL stressCal1
        print*, 'adam'
 1      CONTINUE
        ita = ita + 1
        ita1 = ita1 + 1
        totime = totime + deltat
        !CALL nsMomentumAB
        CALL nsMomentum
        CALL velocityBC
        CALL solidCellBC
        CALL velocityForcing1
        CALL velocityBC
        CALL poissonSolver
        CALL pressureForcing1
        !CALL calculateUcMO
        !CALL calculateUc
        !$acc update host (u, v, w, ut, vt, wt, p, cell)        
        CALL writeOutput  
        !CALL writeOutput1
        !CALL computeSumData 
        !CALL writeComputeSumData         
        !IF(mod(ita1,50000) .eq. 0) THEN
            !CALL computeAvgData
            !CALL writeAvgoutput
        !END IF  
        CALL stressCal2
        !CALL writeStressData
        CALL writeResult        
        IF(ita.lt.itamax) GOTO 1
      END PROGRAM main


         
         
         
         
         
         
         
         
         
         
    
