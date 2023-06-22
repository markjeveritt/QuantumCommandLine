/*
 Code used to produce the data for Figure 10.4.
 
 A third implementation of the method described in.
 
    Quantum-state diffusion with a moving basis: Computing quantum-optical spectra
    Rüdiger Schack, Todd A. Brun, and Ian C. Percival
    Phys. Rev. A 53, 2694 – Published 1 April 1996
    https://doi.org/10.1103/PhysRevA.53.2694
 
 This code took Chapter11_movingBasis_AdaptiveDimension and
 added to it the capability sparse algebraic methods.
 
 Sufficient work was done to reproduce results consistent with results
 in:
 
    Quantum chaos in open systems: a quantum state diffusion analysis
    Todd A. Brun, Ian C. Percival, Rüdiger Schack
    https://doi.org/10.48550/arXiv.quant-ph/9509015
    https://doi.org/10.1088/0305-4470/29/9/020

 However rigorous testing and optimisation has not been done. The code executes
 more slowly than expected and needs profiling to determine why. Maybe for
 real-world use the library should move to using the Accelerate framework
 and parallelisation as well as a better ODE solver.

 The code is not particularly clear and needs to be made clean.
 
 There is much room for improvement in this example.
 */

import Foundation
import Quantum


public func  DuffingQSD_with_moving_basis_sparse() {
    /**
     NB: DuffingOscillatorMovingBasis must be a class as we want to
     pass the reference to it into QuantumStateDiffusionModelMovingBasis
     this is because this deriver will want to print out values from DO
     as it evolves and updates. If we used a struct we would be passing in
     a copy.
    */
//    let DO = DuffingOscillatorMovingBasis(dimension: 40, beta: 0.01)
//    let DO = DuffingOscillatorMovingBasisAdaptiveDimension(beta: 0.01)
    let DO = DuffingOscillatorMovingBasisAdaptiveDimensionSparse(beta: 0.01)
    var QSD = QuantumStateDiffusionModelMovingBasisSparse(forSystem: DO)

    let number_of_periods = 100000
    let steps_per_period = 2000
    let one_period = 2.0 * Double.pi
    let timeIncrement = one_period / Double(steps_per_period)
    
    let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/Duffing/"
    let writeFilename = pathToFile + "PSMB_Sparse_b\(DO.beta).dat"
    
    _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
    let fileHandle = FileHandle(forWritingAtPath: writeFilename)
    
    for _ in 0 ... number_of_periods {
        for _ in 0 ..< steps_per_period {
            QSD.doStep(by: timeIncrement)
        }
        let q = DO.centre.real * (sqrt(2.0) * DO.beta)
        let p = DO.centre.imag * (sqrt(2.0) * DO.beta)
        let outputText = "\(DO.time)\t\(q)\t\(p)\n"
        fileHandle?.write(outputText.data(using: .utf8)!)
    }
}

public typealias SparseMatrixOperator = SparseMatrix<ComplexReal>



// See http://arXiv.org/abs/quant-ph/9506039v2

// Protocol to enable The Dependency Inversion Principle for a general QSD solver.
// See:
//      QuantumStateDiffusionModelMovingBasisSparse

public protocol QSD_SparseTypeWithMovingBasis {
    var time: Real { get set }
    var state: StateVector { get set }
    mutating func displace()
    func getStaticHamiltonian() -> SparseMatrixOperator
    func getDynamicHamiltonian(time: Double) -> SparseMatrixOperator
    func getLindblads() -> [SparseMatrixOperator]
    func getDisplacement() -> ComplexReal
}

// Duffing oscillator struct which solves for QSD
// Not DRY as there should be only one Duffing oscillator struct
// with different solvers for quantum dynamics (and maybe classical too!).
//
// As passes as an argument into QSDTypeWithMovingBasis must be a class (reference type).
public class DuffingOscillatorMovingBasisAdaptiveDimensionSparse: QSD_SparseTypeWithMovingBasis {
    public var state: StateVector
    var spaces: [StateSpace]
    public var time: Real
    
    public let g, beta, Gamma: Double
    public var duffingOscillatorHamiltonian_shifted: [SparseMatrixOperator]
    
    public var lindblad: [SparseMatrixOperator]
    public var lindblad_shifted: [SparseMatrixOperator]
    public var q_shift, p_shift, a_shift: [SparseMatrixOperator]
    public var q, p, a, ad, Id: [SparseMatrixOperator]
    private var displacer: [Displacer]
    
    public var centre: ComplexReal // current origin in phase space
    
    public var currentBasis = 0
    public var basisIncrement = 5
    public var basisTolerance = 1.0e-4

    public init(g: Real = 0.3,
                beta: Real = 0.25,
                Gamma: Real = 0.125) {
        
        self.g = g
        self.beta = beta
        self.Gamma = Gamma
        time = 0.0
        centre = Complex(real: 0, imag: 0)
        
        spaces = []
        a = []       ; ad = []      ; q = []       ; p = [] ; lindblad = []
        a_shift = [] ; q_shift = [] ; p_shift = [] ; lindblad_shifted = []
        duffingOscillatorHamiltonian_shifted = []  ; Id = []
        displacer = []
        for i in 0 ..< 60 {
            let space = StateSpace(dimension: 20 + basisIncrement * i, label: "Space for duffing oscillator")
            displacer.append(Displacer(in: space))
            spaces.append(space)
            // Actually initialise
            let A = space.annihilationOperator
            let Ad = space.creationOperator
            
            a.append(SparseMatrixOperator(from: A))
            ad.append(SparseMatrixOperator(from: Ad))
            
            let Q = (Ad + A)/sqrt(2.0)
            let P = (Ad - A) * Complex(real: 0.0, imag: 1.0) / sqrt(2.0)
            let Lindblad = A * sqrt(2.0 * Gamma)
            q.append(SparseMatrixOperator(from: Q))
            p.append(SparseMatrixOperator(from: P))
            lindblad.append(SparseMatrixOperator(from: Lindblad))
            Id.append(SparseMatrixOperator(from: space.identityOperator))
            
            // To initialise -- will be overwritten
            duffingOscillatorHamiltonian_shifted.append(SparseMatrixOperator(from: space.identityOperator))
            a_shift.append(SparseMatrixOperator(from: A))
            q_shift.append(SparseMatrixOperator(from: Q))
            p_shift.append(SparseMatrixOperator(from: P))
            lindblad_shifted.append(SparseMatrixOperator(from: Lindblad))
        }
        state = spaces[currentBasis].makeNumberState(1)
        
    }
    
    private func sum_squares(array: [ComplexReal] ) -> Real {
        var sum = Real(0)
        
        for element in array {
            sum = sum + element.modulus
        }
        return sum
    }
    
    private func chooseBasis() {
        let dimension = state.elements.count
        let magnitude_last_two_terms = sum_squares(array: Array(state.elements[dimension - 3 ..< dimension]) )
        if ( magnitude_last_two_terms > basisTolerance) { // Grow
            currentBasis = currentBasis + 1
            print("New Basis: \(currentBasis)")

            let zero = ComplexReal(real: 0.0)
            var elements = state.elements
            for _ in 0 ..< basisIncrement {
                elements.append(zero)
            }
            state = spaces[currentBasis].makeVector(from: elements)
            return
        }
        if currentBasis > 1 { // check if we can shrink
            let magnitude_basisIncrement = sum_squares(array: Array(state.elements[dimension - basisIncrement ..< dimension]) )
            if ( magnitude_basisIncrement < basisTolerance * 0.01 ) { // shrink
                currentBasis = currentBasis - 1
                state = spaces[currentBasis].makeVector(from: Array(state.elements[0 ..< dimension - basisIncrement]))
            }
            return
        }
    }
    
    public func displace() {

        chooseBasis() // check what the ode/ito step has done to state
                
        let dAlpha = getDisplacement()
        
        let dQ = sqrt(2) * dAlpha.real
        let dP = sqrt(2) * dAlpha.imag
        let Q = sqrt(2) * centre.real + dQ
        let P = sqrt(2) * centre.imag + dP

        centre = centre + dAlpha
        
        q_shift[currentBasis] = q[currentBasis] + (Id[currentBasis] * Q) //.plus_identity_scaled(by: Q)
        p_shift[currentBasis] = p[currentBasis] + (Id[currentBasis] * P)//.plus_identity_scaled(by: P)
        a_shift[currentBasis] = a[currentBasis] + (Id[currentBasis] * centre)//.plus_identity_scaled(by: centre)
        
        let T = 0.5 * p_shift[currentBasis] * p_shift[currentBasis]
        let qsqr = q_shift[currentBasis] * q_shift[currentBasis]

        let V1 = ( (beta * beta * 0.25) * qsqr * qsqr )
        let V2 = ( -0.5 * qsqr )
        let V3 = 0.5 * Gamma * ( q_shift[currentBasis] * p_shift[currentBasis] + p_shift[currentBasis] * q_shift[currentBasis] )
        
        duffingOscillatorHamiltonian_shifted[currentBasis] = T + V1 + V2 + V3
        
        lindblad_shifted[currentBasis] = lindblad[currentBasis]  + ( Id[currentBasis] * (centre * sqrt(2.0 * Gamma)) ) //.plus_identity_scaled(by: centre * sqrt(2.0 * Gamma) )
        
        let D = displacer[currentBasis].getDispacer(by: -dAlpha) //displacementOperator(alpha: -dAlpha)
        state = D * state
        // fix phase - see eq 16
        let minusIqPminuspQ = Complex(real: 0.0, imag: -( Q * dP - P * dQ) )
        state = state * Complex.exp( minusIqPminuspQ )
        
        chooseBasis() // check what the displacement operator has done to state
    }
    
    func displacementOperator(alpha: ComplexReal) -> MatrixOperator {
        
        let pre_factor = exp( alpha.modulus * alpha.modulus * -0.5 )
        let displacementOperator = Complex(real: pre_factor) *
        
        spaces[currentBasis].exponentialOfScaledCreationOperator(scaleFactor: alpha) *
        spaces[currentBasis].exponentialOfScaledAnnihilationOperator(scaleFactor: -alpha.conjugate )
        
        return displacementOperator
    }
    
    public func getStaticHamiltonian() -> SparseMatrixOperator {
        return duffingOscillatorHamiltonian_shifted[currentBasis]
    }
    
    public func getDynamicHamiltonian(time: Real) -> SparseMatrixOperator {
        return ( ( g / beta ) * cos(time) ) * q_shift[currentBasis]
    }
    
    public func getLindblads() -> [SparseMatrixOperator] {
        return [lindblad_shifted[currentBasis]]
    }
    
    public func getDisplacement() -> ComplexReal {
        // locally this is just the usual `a' as we want changes
        return a[currentBasis].expectationValue(of: state)
    }
}

// In production code would be added to main library.
extension SparseMatrixOperator {
    public func hermitianAdjoint() -> SparseMatrixOperator {
       var output = SparseMatrixOperator(in: self.space)
            for element in self.values {
                output.values.append(CoordinateStorage(value: element.value.conjugate, row: element.col, col: element.row))
            }
        return output
    }
}

// Uses QSD_SparseTypeWithMovingBasis to enable a general QSD solver.
public struct QuantumStateDiffusionModelMovingBasisSparse {
    private var system: QSD_SparseTypeWithMovingBasis
    private let minusi = ComplexReal( real: 0.0, imag: -1.0 )
    
    init(forSystem: QSD_SparseTypeWithMovingBasis) {
        system = forSystem
    }
    
    private func drift(time: Real,
                       psi: StateVector) -> StateVector {
        
        // Schrödinger evolution step
        var psi_out = ( system.getStaticHamiltonian() + system.getDynamicHamiltonian(time: time) ) * (psi * minusi)
        
        // Damping
        let theLinblad = system.getLindblads()
        
        for j in 0 ..< theLinblad.count {
            let lj = theLinblad[j].expectationValue(of: psi)
            psi_out = psi_out + theLinblad[j].hermitianAdjoint() * theLinblad[j] * ( -0.5 * psi )
            psi_out = psi_out + (-0.5 * lj.conjugate * lj) * psi
            psi_out = psi_out + theLinblad[j] * ( lj.conjugate * psi)
        }
        
        return psi_out
    }
    
    private mutating func multiStepIntegrator(by dt: Real) {
        multiStepIvpIntegrator(
            from: system.time,
            to:   system.time + dt,
            first_try_of_stepsize: 0.1 * dt,
            smallest_allowed_value_of_stepsize: 1.0e-10,
            accuracy: 10e-7,
            y: &system.state,
            derivative_function: drift)
        
        system.time += dt
    }
    
    
    public mutating func doStep(by dt: Double) {
        system.displace()
        
        var stochastic_contribution = StateVector(in: system.state.space)
        
        let theLinblad = system.getLindblads()
        for lindbald in theLinblad {
            let psi_1 = lindbald * system.state
            let psi_2 = lindbald.expectationValue(of: system.state) * system.state
            let dEsp = RandomGaussianZeroMeanUnitVariance.nextComplex() * sqrt(dt)
            stochastic_contribution = stochastic_contribution + (psi_1 - psi_2) * dEsp
        }
        
        multiStepIntegrator(by: dt)
        
        let psi = system.state + stochastic_contribution
        // renormalise
        let reciprocal_norm_psi = 1.0 / sqrt( psi.innerProduct(dualVector: psi).real )
        
        system.state = psi * reciprocal_norm_psi
        
        
    }
}

// Percolates parts of the displacement operator to improve performance
// Code needs cleaning

class Displacer {
    private var space: StateSpace
    private var C: [[[Real]]]
    
    init (in space: StateSpace) {
        self.space = space

        let d = space.dimension

        C   = Array(repeating: Array(repeating: Array(repeating: 0.0, count: d), count: d),count: d)
        for n in 0 ..< d  {
            for m in 0 ..< d {
                for k in 0 ... min(m,n) {

                    var temp = 0.5 * ( gammln(n + 1) + gammln(m + 1) )
                    temp -= gammln( n - k + 1) + gammln( m - k + 1) + gammln( k + 1 )

                    C[n][m][k] = exp(temp)
                }
            }
        }
    }
    public func getDispacer(by alpha: ComplexReal) -> MatrixOperator {
        
        var output = space.identityOperator
        let size = space.dimension
        let mzci = -alpha.conjugate

        var iZ   = [ComplexReal]()
        var imZc = [ComplexReal]()
        iZ.append(ComplexReal(real: 1.0))
        imZc.append(ComplexReal(real: 1.0))
        for n in 1 ..< size {
            iZ.append( iZ[n-1] * alpha )   //= Complex.pow(z[i],n);
            imZc.append( imZc[n-1] * mzci)  //= Complex.pow(mzci,n);
        }
        
        let mz    = alpha.modulus
        let tempi = exp(-(mz * mz/(2.0)))
        for n in 0 ..< size  {
           //let cin = C[n]
            for m in 0 ..< size {
              //double cinm[]=cin[m];
               var wkspnm = ComplexReal(real: 0.0)
                for k in 0 ... min(m,n) {
                    wkspnm = wkspnm + ( iZ[n-k] * imZc[m-k] * C[n][m][k] )
                }
                output[n,m] = wkspnm * tempi
            }
        }
        return output
    }
    
    private func gammln(_ x: Int) -> Double {
        return lgamma(Real(x))
    }

}

//
//  Created by M J Everitt on 04/09/2022.
//
