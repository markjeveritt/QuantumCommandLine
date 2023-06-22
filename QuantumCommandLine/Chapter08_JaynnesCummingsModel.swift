/*
 Jaynes Cummings Schrödinger Dynamics driver as described in
 Chapter 8 to produce the atomic inversion seen in Figure 8.2.
 
 Data is produced as tab delimited text.

 See Chapter08_Wigner for code to produce the associated
 Wigner function plot.

 For specifics on model, units and example/test
 parameter values and data see Fig 1 and associated
 text in https://doi.org/10.1103/PhysRevA.79.032328
 also at https://arxiv.org/abs/0710.1983
*/
import Foundation
import Quantum


// MARK: - Set up system
public func  Chapter08_JaynesCummingsModel() {
    
    
    let couplingConstant = 1.0
    let fieldBasisSize = 30
    
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
    
    // MARK: - Set up dynamics
    
    let timeIncrement = 0.01
    let endTime = 35.0
    
    let initialSpinState =  spinSpace.makeVector(from: [1.0, 0.0])
    
    let initialFieldState =  fieldSpace.makeCoherentState(alpha: Complex(modulus: sqrt(15.0), argument: 0.0))
    
    let initalState = totalSpace.tensorProduct(of: initialSpinState, with: initialFieldState)
    
    let jaynseCummingTDSE = TimeIndependentSchrodingerSystem(initialstate: initalState, hamiltonian: jaynesCummingsHamiltonian)
    
    jaynseCummingTDSE.useSparseAlgebra()
    
    let sigmaZextended = SparseMatrix(from: totalSpace.tensorProduct(of: spinSpace.sigmaZ, with: fieldSpace.identityOperator))
    
    
    // MARK: -  now proceed to do calculations.
    var outputText = ""
    
    while jaynseCummingTDSE.time < endTime {
        outputText += "\(jaynseCummingTDSE.time)\t"
        let inversion = sigmaZextended.expectationValue(of: jaynseCummingTDSE.Psi)
        outputText += "\(inversion.real)\n"
        jaynseCummingTDSE.evolve(by: timeIncrement)
    }
    
    do {
        // MARK: - output file setup
        
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/jaynesCummings/"
        let writeFilename = pathToFile + "inversion.dat"
        
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)

        try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)
    } catch {
        errorStream.write("Can not write to output")
    }
    
}
//
//  Created by M J Everitt on 22/01/2022.
//
