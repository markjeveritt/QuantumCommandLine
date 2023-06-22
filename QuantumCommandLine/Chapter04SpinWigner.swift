/*
 In the end this function was not used and is poorly tested. It
 was decided that an in depth discussion of spin Wigner functions
 would have required a level of exposition that went beyond the
 scope of the book.
 
 It is included only to show how a thought process develops as
 it later informed the development of the cover image code.
 
 Please see and use spinWignerParity in UtilityCode.swift over that
 presented here.
 
 Plotting the output can be done in gunplot with the commands
 described in CoverImage.swift
 */

import Foundation
import Quantum

public func Chapter04_SpinWignerFunctionPlot() {

    let spinSpace = StateSpace(dimension: 2, label: "just a Spin Space")
    let psi = spinSpace.makeVector(from: [1.0, 0.0])

    let start_theta = 0.0 //-4.0
    let end_theta   =  Double.pi  // 4.0
    let start_phi = 0.0
    let end_phi   = 2.0 * Double.pi
    
    let Ntheta = 100
    let Nphi = 100
    
    let dTheta = (end_theta - start_theta)/Double(Ntheta)
    let dPhi = (end_phi - start_phi)/Double(Nphi)


    var outputText = ""
    
    // If computing lots of Wigner functions one might want to store and load displacementOperators
    // rather than recomputing every time.
    for thetaCount in 0 ..< Ntheta + 1 {
        let theta = start_theta + Double(thetaCount) * dTheta
        for phiCount in 0 ..< Nphi + 1 {
            let phi = start_phi + Double(phiCount) * dPhi
            let I = Complex(real: 0.0, imag: 1.0)

            let exp_i_simga_z_phi = cos(0.5 * phi) * spinSpace.identityOperator
            + (I * sin(0.5 * phi) ) * spinSpace.sigmaZ
            let exp_i_simga_y_theta = cos(0.5 * theta) * spinSpace.identityOperator
                                    + (I * sin(0.5 * theta) ) * spinSpace.sigmaY
            let spinParity = 0.5 * (spinSpace.identityOperator + sqrt(3.0) * spinSpace.sigmaZ )
            
            let displacementOperator = exp_i_simga_z_phi * exp_i_simga_y_theta
            let Pi = displacementOperator * spinParity * displacementOperator.hermitianAdjoint()
            
            let x = sin(theta)*cos(phi)
            let y = sin(theta)*sin(phi)
            let z = cos(theta)
            let W = Pi.expectationValue(of: psi)

            outputText += "\(x)\t\(y)\t\(z)\t\(W.real)\n"
        }
        outputText += "\n"
    }
    do {
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/WF/"
        let writeFilename = pathToFile + "SpinWigner.dat"
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)

        try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)
    } catch {
        errorStream.write("Can not write to output")
    }
}

//
//  Created by M J Everitt on 16/08/2022.
//
