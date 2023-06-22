/*
 Produces a probability density function for a duffing oscillator with an
 initial condition of a Gaussian distribution in phase space (Figure 2.3).
 
 It uses a Mote-Carlo approach averaging over many trajectories whose initial
 conditions are appropriately sampled.
 
 The code also provide an example of how to use DispatchQueue to perform simple
 parallel computation.
 
 The code contains
 Chapter02_DuffingOscillatorPDF - the top level driver
 DuffingOscillator              - A struct that represents the classical Duffing
                                  oscillator and can evolve its state
 ProbabilityDensityFunctionGrid  - A struct for simplifying the construction of a
                                  probability density function. It contains
                                  a method for outputting the PDF as a tab delimited
                                  String
*/
import Foundation
import Quantum
import GameKit

public func Chapter02_DuffingOscillatorPDF() {

    let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/DuffingOscillator/"
    let writeFilename = pathToFile + "PDF.dat"
    _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)

    // proceed to do calculations.
    
    let number_of_periods = 100
    let number_of_threads = 1000000
    var pdf = ProbabilityDensityFunctionGrid(numberOfXValues: 666,
                                            numberOfYValues: 400,
                                            minimumValueOfXRange: -1.5,
                                            maximumValueOfXRange: 1.5,
                                            minimumValueOfYRange: -0.6,
                                            maximumValueOfYRange: 1.2)


    var setOfOscillators = Array<DuffingOscillator>()
    let random = GKRandomSource()
    let gauisianlimit = 1000000
    let randomGausian = GKGaussianDistribution(randomSource: random, lowestValue: -gauisianlimit, highestValue: gauisianlimit)

    for _ in 0 ..< number_of_threads {
        let initial_position = 1.0 + Double(randomGausian.nextInt())/Double(gauisianlimit * 10)
        let initial_velocity = 1.0 + Double(randomGausian.nextInt())/Double(gauisianlimit * 10)
        setOfOscillators.append(DuffingOscillator(position: initial_position, velocity: initial_velocity))
    }
    for oscilator in setOfOscillators {
        pdf.addPoint(xValue: oscilator.position, yValue: oscilator.velocity)
    }
    
    do {
        try pdf.plotabbleOutput().write(toFile: pathToFile + "pdf_0.dat", atomically: false, encoding: .utf8)
        pdf.resetValuesToZero()
    } catch {
        errorStream.write("Can not write to output")
    }
    
    for period in 1 ... number_of_periods {

        print("\r computing period \(period)")
        
        DispatchQueue.concurrentPerform(iterations: setOfOscillators.count) {
            oscillator in
            setOfOscillators[oscillator].evolveOnePeriod(numberOfSteps: 1000)
        }
        for oscilator in setOfOscillators {
            pdf.addPoint(xValue: oscilator.position, yValue: oscilator.velocity)
        }
        
        
        do {
            try pdf.plotabbleOutput().write(toFile: pathToFile + "pdf_\(period).dat", atomically: false, encoding: .utf8)
            pdf.resetValuesToZero()
        } catch {
            errorStream.write("Can not write to output")
        }
    }
    
}



struct DuffingOscillator {
    // Physical carateristics
    let α, β, δ, γ, ω: Double    // dynamical properties
    var time: Double
    var state: [Double]
    var yerr: [Double]
    // convenience code
    var position: Double { state[0] }
    var velocity: Double { state[1] }
    
    // using standard example values 
    init(position: Double = 1.0,
         velocity: Double = 1.1,
         α: Double = 1.0,
         β: Double = -1.0,
         δ: Double = 0.2,
         γ: Double = 0.3,
         ω: Double = 1.0)
    {
        self.α = α
        self.β = β
        self.δ = δ
        self.γ = γ
        self.ω = ω
        time = 0.0
        state = [position,velocity]
        period = 2.0 * Double.pi / ω
        yerr = [0.0,0.0]
    }

    let period: Double
    public mutating func evolveOnePeriod(numberOfSteps: Int) {
// TODO: Why does this not work? Is it just non-linear dynmaics on an integrator issue?
//        multiStepIvpIntegrator(from: time,
//                               to: time + period,
//                               first_try_of_stepsize: 0.01,
//                               smallest_allowed_value_of_stepsize: period / 10000.0,
//                               accuracy: 1.0e-8,
//                               y: &state,
//                               derivative_function: equationOfMotion)
//        time = time + period
        let dt = period/Double(numberOfSteps)
        for _ in 0 ..< numberOfSteps {
            evolve(by: dt)
        }
    }
    
    private func equationOfMotion(t: Double, y: [Double]) -> [Double] {
        return  [ velocity ,
                  γ * cos (ω * t) - δ * velocity - position * (β + α * position * position ) ]
    }
    
    public mutating func evolve(by dt: Double) {
        
//  state = doRungeKuttaStep(t: time,
//                                 h: dt,
//                                 y: state,
//                                 derivative_function: equationOfMotion)
//  let dydx = equationOfMotion(t: time, y: state)
        
//  fifthOrderCashKarpRungeKutta(t: time, h: dt, y: state, yout: &state, yerr: &yerr, derivative_function: equationOfMotion(t:y:), dydx: dydx)
        
        multiStepIvpIntegrator(from: time,
                               to: time + dt,
                               first_try_of_stepsize: dt/10.0,
                               smallest_allowed_value_of_stepsize: 1.0e-8,
                               accuracy: 1.0e-4,
                               y: &state,
                               derivative_function: equationOfMotion)
        time += dt
    }
}

struct ProbabilityDensityFunctionGrid {
    
    private let numberOfXValues: Int
    private let minX: Double
    private let maxX: Double
    private let minY: Double
    private let maxY: Double
    private let numberOfYValues: Int
    private var grid: [Int]
    private let xticks: [Double]
    private let yticks: [Double]
    private var normalisation = 0
    
    public init (numberOfXValues nx: Int,
                 numberOfYValues ny: Int,
                 minimumValueOfXRange: Double,
                 maximumValueOfXRange: Double,
                 minimumValueOfYRange: Double,
                 maximumValueOfYRange: Double
    ) {
        numberOfXValues = nx
        numberOfYValues = ny
        minX = minimumValueOfXRange
        maxX = maximumValueOfXRange
        minY = minimumValueOfYRange
        maxY = maximumValueOfYRange
        
        grid = Array.init(repeating: 0, count: numberOfXValues * numberOfYValues)
        
        var xticks_temp = Array.init(repeating: 0.0, count: numberOfXValues)
        let dx = (maxX - minX) / Double(nx)
        for i in 0 ..< numberOfXValues {
            xticks_temp[i] = minimumValueOfXRange + Double(i) * dx
        } // doing this as never want to change the tick values once set
        xticks = xticks_temp

        var yticks_temp = Array.init(repeating: 0.0, count: numberOfYValues)
        let dy = (maxY - minY) / Double(ny)
        for i in 0 ..< numberOfYValues {
            yticks_temp[i] = minimumValueOfYRange + Double(i) * dy
        }
        yticks = yticks_temp
    }
    
    public mutating func resetValuesToZero() {
        grid = grid.map( { _ in 0 } )
        normalisation = 0
    }
    
    public mutating func addPoint(xValue x: Double, yValue y: Double) {
        
        let xIndex = xticks.lastIndex(where: {$0 <= x} ) ?? -1
        let yIndex = yticks.lastIndex(where: {$0 <= y} ) ?? -1
        
        if (xIndex >= 0) && (yIndex >= 0) {
            self[xIndex,yIndex] += 1
            normalisation += 1
        } else {
            print("Warining: point outside PDF")
        }
        
    }
    
    subscript(row: Int, col: Int) -> Int {
        get {
            let index = atIndex(row: row, column: col)
            return grid[index]
        }
        set {
            let index = atIndex(row: row, column: col)
            grid[index] = newValue
        }
    }
    
    func atIndex(row: Int, column: Int) -> Int {
        return (row * numberOfYValues) + column
    }
    
    public func plotabbleOutput() -> String {
        var output = ""
        
        for x in 0 ..< xticks.count {
            for y in 0 ..< yticks.count {
                output.append("\(xticks[x])\t\(yticks[y])\t\(Double(self[x,y]) / Double(normalisation))")
                output.append("\n")
            }
            output.append("\n")
        }
        return output
    }
}

//
//  Created by M J Everitt on 22/01/2022.
//
