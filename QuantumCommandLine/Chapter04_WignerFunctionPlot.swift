/*
 Code used to make the Wigner functions in Figures 4.1 and 4.2. over the
 phase space of position (q) and momentum (p).
 
 Note that the method implement here of using the displaced parity operator
 is far from optimised for performance (we dont even account for the fact that
 the parity operator is diagonal - although improvements in the way the
 library could take that into account if e.g. operators had an isDiagonal
 attribute)
 */

import Foundation
import Quantum

public func Chapter04_WignerFunctionPlot() {

    let stateSpace = StateSpace(dimension: 40, label: "just a stateSpace")
    //let psi = stateSpace.makeCoherentState( alpha: Complex(real: 0.0) )
    let psi = stateSpace.makeSqueezedVacuumState( r: 0.5, phi: Double.pi/2.0 )
    //    let psi = stateSpace.makeCoherentState( alpha: Complex(real: 2.0) )
    //    let psi = stateSpace.makeCoherentState( alpha: Complex(real: 2.0) ) + stateSpace.makeCoherentState(alpha: Complex(real: -2.0))

    let start_q = -2.0 //-4.0
    let end_q   =  2.0 // 4.0
    let start_p = -2.0
    let end_p   =  2.0
    
    let Nq = 100
    let Np = 100
    
    let dq = (end_q - start_q)/Double(Nq)
    let dp = (end_p - start_p)/Double(Np)


    var outputText = ""
    
    // If computing lots of Wigner functions one might want to store and load displacementOperators
    // rather than recomputing every time.
    for qCount in 0 ..< Nq {
        let q = start_q + Double(qCount) * dq
        for pCount in 0 ..< Np {
            let p = start_p + Double(pCount) * dp
            let twoalpha = Complex(real: q, imag: p)
            let pre_factor = exp( twoalpha.modulus * twoalpha.modulus * -0.5 )
            // https://journals.aps.org/pra/abstract/10.1103/PhysRevA.50.4488 D(2a)Pi = D(a) Pi D^dag(a)
            let displacementOperator = Complex(real: pre_factor) *
                    stateSpace.exponentialOfScaledCreationOperator(scaleFactor: twoalpha) *
                    stateSpace.exponentialOfScaledAnnihilationOperator(scaleFactor: -(twoalpha.conjugate) )
            let Pi = displacementOperator * stateSpace.parityOperator
            
            let W = Pi.expectationValue(of: psi)
            outputText += "\(q)\t\(p)\t\(W.real)\n"
        }
        outputText += "\n"
    }
    do {
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/WF/"
        let writeFilename = pathToFile + "Wigner.dat"
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)

        try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)
    } catch {
        errorStream.write("Can not write to output")
    }
}


//
//  Created by M J Everitt on 08/08/2022.
//
