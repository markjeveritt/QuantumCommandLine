/*
 Copies code from Chapter08_JaynesCummingsModel.swift and
 adds the functionality to print out the Wigner function
 in Figure 8.2.
 
 This is an example of bad practice that is often seen in
 scientific coding where working code for solving one problem
 is used as a template for solving a related problem. The issue
 is that this is not DRY as the process results in repetition of
 quite a few lines of code.
 
 This could have been avoided, or at least much reduced, by
 encapsulating the code in Chapter08_JaynesCummingsModel in
 sensible structs and funcs (such as a JaynesCummingsModel
 struct). This is left as an exercise.
*/

import Foundation
import Quantum

public func  chapter08Wigner() {
    // MARK: -  Code copied from Chapter08_JaynesCummingsModel from here
    let couplingConstant = 1.0
    let fieldBasisSize = 40
    
    let fieldSpace = StateSpace(dimension: fieldBasisSize, label: "Field")
    
    let fieldHamiltonian = fieldSpace.numberOperator + ( 0.5 * fieldSpace.identityOperator )
    
    let spinSpace  = StateSpace(dimension: 2, label: "Atom")
    let spinHamiltonian = spinSpace.sigmaZ
    
    let totalSpace = StateSpace(tensorProductOf: fieldSpace, spinSpace, label: "JC System")
    
    let interactionTerm = couplingConstant * (
        totalSpace.tensorProduct(
            of: fieldSpace.creationOperator,
            with: spinSpace.sigmaMinus
        )
        +
        totalSpace.tensorProduct(
            of: fieldSpace.annihilationOperator,
            with: spinSpace.sigmaPlus
        )
    )
    
    let fieldHamiltonianExtension = totalSpace.tensorProduct(of: fieldHamiltonian, with: spinSpace.identityOperator)
    let spinHamiltonianExtension = totalSpace.tensorProduct(of: spinHamiltonian, with: fieldSpace.identityOperator)
    
    let jaynesCummingsHamiltonian = fieldHamiltonianExtension + spinHamiltonianExtension + interactionTerm
    
    
    let timeIncrement = 0.01
    // MARK: - changed end time to be in mid cat
    let endTime = 10.0
    
    let initialSpinState =  spinSpace.makeVector(from: [1.0, 0.0])
    
    let initialFieldState =  fieldSpace.makeCoherentState(alpha: Complex(modulus: sqrt(15.0), argument: 0.0))
    
    let initalState = totalSpace.tensorProduct(of: initialSpinState, with: initialFieldState)
    
    let jaynseCummingTDSE = TimeIndependentSchrodingerSystem(initialstate: initalState, hamiltonian: jaynesCummingsHamiltonian)
    
    jaynseCummingTDSE.useSparseAlgebra()
    
//    let sigmaZextended = SparseMatrix(from: totalSpace.tensorProduct(of: spinSpace.sigmaZ, with: fieldSpace.identityOperator))
    
    
    // MARK: -  now proceed to do calculations.
    
    while jaynseCummingTDSE.time < endTime {
//       // let inversion = sigmaZextended.expectationValue(of: jaynseCummingTDSE.Psi)
        jaynseCummingTDSE.evolve(by: timeIncrement)
    }


    // MARK: -  to Here

    // MARK: -  now proceed to do compute Wigner function.

    let start_q = -5.0 //-4.0
    let end_q   =  3.0 // 4.0
    let start_p = -4.0
    let end_p   =  5.0
    
    let Nq = 800
    let Np = 800
    
    let dq = (end_q - start_q)/Double(Nq)
    let dp = (end_p - start_p)/Double(Np)


    var outputText = ""
    
    // If computing lots of Wigner functions one might want to store and load displacementOperators
    // rather than recomputing every time.
    for qCount in 0 ..< Nq {
        let q = start_q + Double(qCount) * dq
        for pCount in 0 ..< Np {
            let p = start_p + Double(pCount) * dp
            let twoalpha = 2.0 * Complex(real: q, imag: p)
            let pre_factor = exp( twoalpha.modulus * twoalpha.modulus * -0.5 )
            // https://journals.aps.org/pra/abstract/10.1103/PhysRevA.50.4488 D(2a)Pi = D(a) Pi D^dag(a)
            let displacementOperator = Complex(real: pre_factor) *
            fieldSpace.exponentialOfScaledCreationOperator(scaleFactor: twoalpha) *
            fieldSpace.exponentialOfScaledAnnihilationOperator(scaleFactor: -(twoalpha.conjugate) )
            let Pi = displacementOperator * fieldSpace.parityOperator
            
            let PiExtended = totalSpace.tensorProduct(of: Pi, with: spinSpace.identityOperator)
            let W = PiExtended.expectationValue(of: jaynseCummingTDSE.Psi )
            outputText += "\(q)\t\(p)\t\(W.real)\n"
        }
        outputText += "\n"
    }
    do {
        // MARK: - output file setup
        
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/Cover/"
        let writeFilename = pathToFile + "WignerFunctionJC.dat"
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)

        try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)

    } catch {
        errorStream.write("Can not write to output")
    }
}


//
//  Created by M J Everitt on 15/08/2022.
//
