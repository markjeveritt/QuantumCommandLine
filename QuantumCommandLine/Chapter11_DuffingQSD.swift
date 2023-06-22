/*
 
 This was test code that was never used to produce any of figures in the book.
 
 The main issue was that it became rapidly apparent that even on modern
 computers this algorithm was insufficiently performant to be of use in
 achieving a good classical limit.
 
 Instead of alternative approaches were adopted.
 
 Three other versions of the code are included as an example of testing out an idea
 and developing a solution in this order:
    - Chapter11_QSD_Moving_Basis
    - Chapter11_movingBasis_AdaptiveDimension
    - Chapter_11_DuffingQSD_with_moving_basis_sparse
 */

import Foundation
import Quantum
public func  DuffingQSD() {
    
    let duffingOscilatorBasisSize = 180
    print("Setting up for basis: \(duffingOscilatorBasisSize)")

    let duffingOscillatorSpace = StateSpace(dimension: duffingOscilatorBasisSize, label: "Duffing")
    let I = Complex(real: 0.0, imag: 1.0)
    
    let a = duffingOscillatorSpace.annihilationOperator
    let a_dagger = duffingOscillatorSpace.creationOperator
    
    let q = (a_dagger + a)/sqrt(2.0)
    let p = (a_dagger - a) * I / sqrt(2.0)

    let initialState =  duffingOscillatorSpace.makeCoherentState(alpha: Complex(modulus: sqrt(3.0), argument: 0.0))
    // parameters taken from Iain Perceval's book quantum state diffusion p92
    let g = 0.3
    let beta = 1.0
    let Gamma = 0.125
    let T = 0.5 * p * p
    let qsqr = q * q
    let V1 = ( (beta * beta * 0.25) * qsqr * qsqr )
    let V2 = ( -0.5 * qsqr )
    let V3 = (0.5 * Gamma) * (q * p + p * q)
    let duffingOscillatorHamiltonian = T + V1 + V2 + V3

    func drive(_ t: Double ) -> Double {
        return ( g / beta ) * cos(t)
    }
    print("making QSD runner")

    let QSD = QuantumStateDiffusionModelSparseAlgebra(
                    initialState: initialState,
                    hamiltonian: duffingOscillatorHamiltonian,
                    driveFunction: drive,
                    driveHamiltonianOperator: q,
                    lindblads: (2.0 * sqrt(Gamma)) * a
                    )
        
    // MARK: - Set up dynamics

    let number_of_periods = 10000
    let steps_per_period = 800
    let timeIncrement = 2.0 * Double.pi / Double(steps_per_period)
    print("Dynamics steps_per_period \(steps_per_period) and \(number_of_periods) periods")
    // MARK: -  now proceed to do dynamics.
    var outputText = ""
    do {
        // MARK: - output file setup
        let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/Duffing/"
        let writeFilename = pathToFile + "PS_b1.0.dat"
        
        _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
        let fileHandle = FileHandle(forWritingAtPath: writeFilename)
        
        
        for i in 0 ... number_of_periods {
            for _ in 0 ..< steps_per_period {
                QSD.doStep(by: timeIncrement)
            }
            let aValue = a.expectationValue(of: QSD.psi)
            print("done \(i) of \(number_of_periods)")
            outputText = "\(aValue.real)\t\(aValue.imag)\n"
            //try outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8)
            try fileHandle?.write(outputText.data(using: .utf8)!)
            //(contentsOf: outputText.data(using: .utf8) ?? 0 )
        }
        fileHandle!.closeFile()
    } catch {
        errorStream.write("Can not write to output")
    }
}
// basic implementation of a quantum jumps model.
// https://arxiv.org/abs/quant-ph/9701024

open class QuantumStateDiffusionModelSparseAlgebra {
    public var psi: StateVector
    private var minus_i_H: SparseMatrix<ComplexReal>
    public var time: Real
    private var theLinblad: [SparseMatrix<ComplexReal>]
    private var LdagL: [SparseMatrix<ComplexReal>]
    private var driveFunction: (Real) -> Real
    private var minus_i_H_time: SparseMatrix<ComplexReal>
    private let minusi = ComplexReal( real: 0.0, imag: -1.0 )
    
    public init(
                initialState: StateVector,
                hamiltonian: MatrixOperator,
                driveFunction: @escaping (Real) -> Real ,
                driveHamiltonianOperator: MatrixOperator,
                lindblads: MatrixOperator...
    ) {
        time = 0.0
        psi = initialState
        minus_i_H = SparseMatrix(from: minusi * hamiltonian)
        theLinblad = []
        LdagL = []
        for lindblad in lindblads {
            theLinblad.append(SparseMatrix(from: lindblad))

            let Ldag_L = lindblad.hermitianAdjoint() * lindblad
            LdagL.append(SparseMatrix(from: Ldag_L))
        }
        self.driveFunction = driveFunction
        minus_i_H_time = SparseMatrix(from: minusi * driveHamiltonianOperator)
    }

    private func drift(
        time: Real,
        psi: StateVector) -> StateVector {
        
        // Schrödinger evolution step
        var psi_out = ( minus_i_H + minus_i_H_time * driveFunction(time) ) * psi
        
        // Drift term
        for j in 0 ..< theLinblad.count {
            let lj = theLinblad[j].expectationValue(of: psi)
            psi_out = psi_out + LdagL[j] * ( -0.5 * psi )
            psi_out = psi_out + (-0.5 * lj.conjugate * lj) * psi
            psi_out = psi_out + theLinblad[j] * ( lj.conjugate * psi)
        }
        return psi_out
    }
        
    private func multiStepIntegrator(by dt: Real) {
        multiStepIvpIntegrator(
            from: time,
            to:   time + dt,
            first_try_of_stepsize: dt,
            smallest_allowed_value_of_stepsize: 1.0e-8,
            accuracy: 10e-6,
            y: &psi,
            derivative_function: drift)
        time += dt
    }
    
    
    public func doStep(by dt: Double) {
        
        var stochastic_contribution = StateVector(in: psi.space)
        
        for lindbald in theLinblad {
            let psi_1 = lindbald * psi
            let psi_2 = lindbald.expectationValue(of: psi) * psi
            let dEsp = RandomGaussianZeroMeanUnitVariance.nextComplex() * sqrt(dt)
            stochastic_contribution = stochastic_contribution + (psi_1 - psi_2) * dEsp
        }
        
        multiStepIntegrator(by: dt)

        psi = psi + stochastic_contribution

        // renormalise
        let reciprocal_norm_psi = 1.0 / sqrt( psi.innerProduct(dualVector: psi).real )
        psi = psi * reciprocal_norm_psi
        
    }
}


//
//  Created by M J Everitt on 22/08/2022.
//

