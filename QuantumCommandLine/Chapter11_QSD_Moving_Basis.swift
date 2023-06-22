/*
 A first implementation of the method described in.
 
    Quantum-state diffusion with a moving basis: Computing quantum-optical spectra
    Rüdiger Schack, Todd A. Brun, and Ian C. Percival
    Phys. Rev. A 53, 2694 – Published 1 April 1996
    https://doi.org/10.1103/PhysRevA.53.2694
 
 Code ran too slowly to be useful and was superseded by
 Chapter11_movingBasis_AdaptiveDimension.swift
 */

import Foundation
import Quantum

public func  DuffingQSD_with_moving_basis() {
    /**
     NB: DuffingOscillatorMovingBasis must be a class as we want to
     pass the reference to it into QuantumStateDiffusionModelMovingBasis
     this is because this deriver will want to print out values from DO
     as it evolves and updates. If we used a struct we would be passing in
     a copy.
    */
    let DO = DuffingOscillatorMovingBasis(dimension: 40, beta: 0.25)
    var QSD = QuantumStateDiffusionModelMovingBasis(forSystem: DO)

    let number_of_periods = 100000
    let steps_per_period = 2000
    let one_period = 2.0 * Double.pi
    let timeIncrement = one_period / Double(steps_per_period)
    
    let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/Duffing/"
    let writeFilename = pathToFile + "PSMB_Simple_b\(DO.beta).dat"
    
    _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
    let fileHandle = FileHandle(forWritingAtPath: writeFilename)
    
    for _ in 0 ... number_of_periods {
        
        for _ in 0 ..< steps_per_period {
            QSD.doStep(by: timeIncrement)
            let q = DO.centre.real * (sqrt(2.0) * DO.beta)
            let p = DO.centre.imag * (sqrt(2.0) * DO.beta)
            //print("\(DO.time * 0.5 / Double.pi)")
            let outputText = "\(DO.time)\t\(q)\t\(p)\n"
            fileHandle?.write(outputText.data(using: .utf8)!)
        }
        
    }
}

// http://arXiv.org/abs/quant-ph/9506039v2
// As passed as an argument into QSDTypeWithMovingBasis must be a class (reference type).
public class DuffingOscillatorMovingBasis: QSDTypeWithMovingBasis {
    var state: StateVector
    var space: StateSpace
    var time: Real
    
    public let g, beta, Gamma: Double
    public var duffingOscillatorHamiltonian_shifted: MatrixOperator
    
    public var lindblad: MatrixOperator
    public var lindblad_shifted: MatrixOperator
    public var q_shift, p_shift, a_shift: MatrixOperator
    public let q, p, a, ad: MatrixOperator
    
    public var centre: ComplexReal // current origin in phase space
    
    public init(dimension: Int,
                g: Real = 0.3,
                beta: Real = 0.25,
                Gamma: Real = 0.125) {
        
        self.g = g
        self.beta = beta
        self.Gamma = Gamma
        time = 0.0
        centre = Complex(real: 0, imag: 0)
        space = StateSpace(dimension: dimension, label: "Space for duffing oscillator")
        state = space.makeNumberState(0)
        
        // Actually initialise
        a = space.annihilationOperator
        ad = space.creationOperator
        q = (ad + a)/sqrt(2.0)
        p = (ad - a) * Complex(real: 0.0, imag: 1.0) / sqrt(2.0)
        lindblad = space.annihilationOperator * sqrt(2.0 * Gamma)
        
        // To initialise -- will be overwritten
        duffingOscillatorHamiltonian_shifted = space.identityOperator
        a_shift = a
        q_shift = q
        p_shift = p
        lindblad_shifted = lindblad
    }
    
    public func displace() {
        let dAlpha = getDisplacement()
        
        let dQ = sqrt(2) * dAlpha.real
        let dP = sqrt(2) * dAlpha.imag
        let Q = sqrt(2) * centre.real + dQ
        let P = sqrt(2) * centre.imag + dP

        centre = centre + dAlpha
        
        q_shift = q.plus_identity_scaled(by: Q)
        p_shift = p.plus_identity_scaled(by: P)
        a_shift = a.plus_identity_scaled(by: centre)
        
        let T = 0.5 * p_shift * p_shift
        let qsqr = q_shift * q_shift

        let V1 = ( (beta * beta * 0.25) * qsqr * qsqr )
        let V2 = ( -0.5 * qsqr )
        let V3 = 0.5 * Gamma * ( q_shift * p_shift + p_shift * q_shift )
        
        duffingOscillatorHamiltonian_shifted = T + V1 + V2 + V3
        
        lindblad_shifted = lindblad.plus_identity_scaled(by: centre * sqrt(2.0 * Gamma) )
        
        let D = displacementOperator(alpha: -dAlpha)
        state = D * state
        // fix phase - see eq 16
        let minusIqPminuspQ = Complex(real: 0.0, imag: ( Q * dP - P * dQ) )
        state = state * Complex.exp( minusIqPminuspQ )
        
    }
    
    func displacementOperator(alpha: ComplexReal) -> MatrixOperator {
        
        let pre_factor = exp( alpha.modulus * alpha.modulus * -0.5 )
        let displacementOperator = Complex(real: pre_factor) *
        
        space.exponentialOfScaledCreationOperator(scaleFactor: alpha) *
        space.exponentialOfScaledAnnihilationOperator(scaleFactor: -alpha.conjugate )
        
        return displacementOperator
    }
    
    func getStaticHamiltonian() -> MatrixOperator {
        return duffingOscillatorHamiltonian_shifted
    }
    
    func getDynamicHamiltonian(time: Real) -> MatrixOperator {
        return ( ( g / beta ) * cos(time) ) * q_shift
    }
    
    func getLindblads() -> [MatrixOperator] {
        return [lindblad_shifted]
    }
    
    func getDisplacement() -> ComplexReal {
        // locally this is just the usual `a' as we want changes
        return a.expectationValue(of: state)
    }
}


protocol QSDTypeWithMovingBasis {
    var time: Real { get set }
    var state: StateVector { get set }
    mutating func displace()
    func getStaticHamiltonian() -> MatrixOperator
    func getDynamicHamiltonian(time: Double) -> MatrixOperator
    func getLindblads() -> [MatrixOperator]
    func getDisplacement() -> ComplexReal
}
public struct QuantumStateDiffusionModelMovingBasis {
    //public var psi: StateVector
    
    private var system: QSDTypeWithMovingBasis
    
    private let minusi = ComplexReal( real: 0.0, imag: -1.0 )
    
    init(forSystem: QSDTypeWithMovingBasis) {
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
        
         
        // Failed for beta = 0.01
        // Unable to determine why in timescale available
        // May just need a better integration routine or more basis states than
        // computationally acceptable.
         
        multiStepIvpIntegrator(
            from: system.time,
            to:   system.time + dt,
            first_try_of_stepsize: 0.1 * dt,
            smallest_allowed_value_of_stepsize: 1.0e-10,
            accuracy: 10e-7,
            y: &system.state,
            derivative_function: drift)
        
            // First attempt used a simple integrator.
            // system.state = doRungeKuttaStep(t: system.time,
            //                                 h: dt,
            //                                 y: system.state,
            //                                 derivative_function: drift)
        
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

//
//  Created by M J Everitt on 27/08/2022.
//
