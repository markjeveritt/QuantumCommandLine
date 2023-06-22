/*
 A second implementation of the method described in.
 
    Quantum-state diffusion with a moving basis: Computing quantum-optical spectra
    Rüdiger Schack, Todd A. Brun, and Ian C. Percival
    Phys. Rev. A 53, 2694 – Published 1 April 1996
    https://doi.org/10.1103/PhysRevA.53.2694
 
 This code took Chapter11_QSD_Moving_Basis.swift and
 added to it the capability to changes the number of biases
 states dynamically within the simulation.
 
 Code ran too slowly to be useful and was superseded by
 Chapter_11_DuffingQSD_with_moving_basis_sparse.swift
 */

import Foundation
import Quantum

// As passed as an argument into QSDTypeWithMovingBasis must be a class (reference type).

public class DuffingOscillatorMovingBasisAdaptiveDimension: QSDTypeWithMovingBasis {
    var state: StateVector
    var spaces: [StateSpace]
    var time: Real
    
    public let g, beta, Gamma: Double
    public var duffingOscillatorHamiltonian_shifted: [MatrixOperator]
    
    public var lindblad: [MatrixOperator]
    public var lindblad_shifted: [MatrixOperator]
    public var q_shift, p_shift, a_shift: [MatrixOperator]
    public var q, p, a, ad: [MatrixOperator]
    
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
        duffingOscillatorHamiltonian_shifted = []
        
        for i in 0 ..< 40 {
            let space = StateSpace(dimension: 20 + basisIncrement * i, label: "Space for duffing oscillator")
            spaces.append(space)
            // Actually initialise
            let A = space.annihilationOperator
            let Ad = space.creationOperator
            
            a.append(A)
            ad.append(Ad)
            
            let Q = (Ad + A)/sqrt(2.0)
            let P = (Ad - A) * Complex(real: 0.0, imag: 1.0) / sqrt(2.0)
            let Lindblad = A * sqrt(2.0 * Gamma)
            q.append(Q)
            p.append(P)
            lindblad.append(Lindblad)
            
            // To initialise -- will be overwritten
            duffingOscillatorHamiltonian_shifted.append(space.identityOperator)
            a_shift.append(A)
            q_shift.append(Q)
            p_shift.append(P)
            lindblad_shifted.append(Lindblad)
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
        
        q_shift[currentBasis] = q[currentBasis].plus_identity_scaled(by: Q)
        p_shift[currentBasis] = p[currentBasis].plus_identity_scaled(by: P)
        a_shift[currentBasis] = a[currentBasis].plus_identity_scaled(by: centre)
        
        let T = 0.5 * p_shift[currentBasis] * p_shift[currentBasis]
        let qsqr = q_shift[currentBasis] * q_shift[currentBasis]

        let V1 = ( (beta * beta * 0.25) * qsqr * qsqr )
        let V2 = ( -0.5 * qsqr )
        let V3 = 0.5 * Gamma * ( q_shift[currentBasis] * p_shift[currentBasis] + p_shift[currentBasis] * q_shift[currentBasis] )
        
        duffingOscillatorHamiltonian_shifted[currentBasis] = T + V1 + V2 + V3
        
        lindblad_shifted[currentBasis] = lindblad[currentBasis].plus_identity_scaled(by: centre * sqrt(2.0 * Gamma) )
        
        let D = displacementOperator(alpha: -dAlpha)
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
    
    func getStaticHamiltonian() -> MatrixOperator {
        return duffingOscillatorHamiltonian_shifted[currentBasis]
    }
    
    func getDynamicHamiltonian(time: Real) -> MatrixOperator {
        return ( ( g / beta ) * cos(time) ) * q_shift[currentBasis]
    }
    
    func getLindblads() -> [MatrixOperator] {
        return [lindblad_shifted[currentBasis]]
    }
    
    func getDisplacement() -> ComplexReal {
        // locally this is just the usual `a' as we want changes
        return a[currentBasis].expectationValue(of: state)
    }
}

//
//  Created by M J Everitt on 27/08/2022.
//
