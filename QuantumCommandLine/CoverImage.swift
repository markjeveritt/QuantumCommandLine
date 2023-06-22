/*
 Generates the data to produce the cover image as described in
 section 4.3.6. Generalising Wigner Functions to Spins.
 
 This contains only the main driver code and the function on
 which it relies can be found in UtilityCode.swift
 
 The plots was made using gunplot with pm3d and lighting. Some
 relevant gnuplot commands for producing simple plots from the
 output of this driver are:
 
    set pm3d lighting primary 0.5 specular 0.2 spec2 0
    set pm3d at s
    unset sur
    set pm3d depthorder
    set pm3d lighting primary 0.1
    set palette defined ( 0.2 "dark-blue", 0.5 "grey90", 1 "dark-red" )
    set view equal xyz
    splot "WignerFunctionSpinCat.dat" u 1:2:3:4
 
 */

import Foundation
import Quantum

public func  coverImage() {
    
    let fieldBasisSize = 20
    let fieldSpace = StateSpace(dimension: fieldBasisSize, label: "Field")
    let spinSpace  = StateSpace(dimension: 2, label: "Atom")
    let totalSpace = StateSpace(tensorProductOf: fieldSpace, spinSpace, label: "JC System")

    let initialSpinState =  spinSpace.makeVector(from: [1.0, 0.0])
    let initialFieldState =  fieldSpace.makeCoherentState(alpha: Complex(modulus: sqrt(15.0), argument: 0.0))
    let a = 1.0/sqrt(2.0)
    let spinUp   = spinSpace.makeVector(from: [1, 0])
    let spinDown = spinSpace.makeVector(from: [0, 1])

    let alpha2_5      = fieldSpace.makeCoherentState(alpha: Complex(modulus: 2, argument: 0.0))
    let alphaMinus2_5 = fieldSpace.makeCoherentState(alpha: Complex(modulus: -2, argument: 0.0))

    let psi = ( totalSpace.tensorProduct(of: spinUp, with: alpha2_5) +
                totalSpace.tensorProduct(of: spinDown, with: alphaMinus2_5) ) / sqrt(2.0)
 //   let psi = ( totalSpace.tensorProduct(of: spinUp + spinDown, with: alpha2_5 + alphaMinus2_5) ) / 2.0

    let start_q = -3.5, end_q = 3.5, start_p = -1.5, end_p = 1.5
    let start_theta = 0.0, end_theta = Double.pi, start_phi = 0.0, end_phi = 2.0 * Double.pi

    let Nq = 80
    let Ntheta = 80, Nphi = 80

    
    let dq = (end_q - start_q) / Double(Nq)
    let R = dq / 2.0
    let dp = dq //(end_p - start_p)/Double(Np)
    let Np = Int( (end_p - start_p) / dq )
    let dTheta = (end_theta - start_theta)/Double(Ntheta)
    let dPhi = (end_phi - start_phi)/Double(Nphi)

    
    var outputText = ""
    
    for qCount in 0 ... Nq {
        let q = start_q + Double(qCount) * dq
        for pCount in 0 ... Np {
            
            let p = start_p + Double(pCount) * dp
            let PiField = fieldWignerParity(q: q,
                                            p: p,
                                            fieldSpace: fieldSpace)
            let PiFieldExtended = totalSpace.tensorProduct(of:   PiField,
                                                           with: spinSpace.identityOperator)
            let W = (PiFieldExtended.expectationValue(of: psi)).real
            
            for thetaCount in 0 ... Ntheta  {
                let theta = start_theta + Double(thetaCount) * dTheta
                for phiCount in 0 ... Nphi  {

                    let phi = start_phi + Double(phiCount) * dPhi
                    
                    let piSpin = spinWignerParity(theta: theta,
                                                  phi: phi,
                                                  spinSpace: spinSpace)
                    
                    let PiSpinTotal = totalSpace.tensorProduct(of:   piSpin,
                                                               with: PiField)
                    
                    let x = q + R * sin(theta) * cos(phi)
                    let y = p + R * sin(theta) * sin(phi)
                    let z = 1.2 * W + R * cos(theta)
                    
                    let Wspin = PiSpinTotal.expectationValue(of: psi )
                    outputText += "\(x)\t\(y)\t\(z)\t\(Wspin.real)\n"
                }
                outputText += "\n"
            }
            outputText += "\n"
        }
        outputText += "\n"
    }
    do {
        // MARK: - output file setup
        
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/Cover/"
        let writeFilename = pathToFile + "WignerFunctionSpinCat.dat"
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
        
        try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)
        
    } catch {
        errorStream.write("Can not write to output")
    }
}

//
//  Created by M J Everitt on 15/08/2022.
//
