#!/usr/bin/env ruby
# Copyright 2023, Christopher Dilks
# Subject to the terms in the LICENSE file found in the top-level directory.

require 'optparse'
require 'ostruct'
require 'fileutils'

###################################
# constants:
EtaTestValues = {
  :ideal => 2.0,
  :min   => 1.6,
  :max   => 3.5,
}
IdealEnergy   = 12.0 # [GeV]
IdealParticle = 'pi+'
###################################

# setup
# ---------------------------------------------------

# default opt
opt = OpenStruct.new
opt.sim_mode      = ''
opt.sim_file      = 'out/sim.edm4hep.root'
opt.rec_file      = 'out/rec.edm4hep.root'
opt.ana_file      = 'out/ana.edm4hep.root'
opt.run_sim       = true
opt.run_rec       = true
opt.run_ana       = true
opt.run_rec_down  = false
opt.run_ana_down  = false
opt.num_events    = 10
opt.benchmark_exe = 'benchmark_rich_reconstruction'
opt.algos         = Array.new
opt.plots         = Array.new
opt.verbosity     = 0
opt.dry_run       = false
opt.using_ci      = false

# available simulation modes
avail_sim_modes = [
  'fixedEtaIdeal',
  'fixedEtaMin',
  'fixedEtaMax',
]

# parse options
required_set = false
OptionParser.new do |o|
  o.banner = "USAGE: #{$0} [OPTIONS]..."
  o.separator ''
  o.separator 'required options, one of either:'.upcase
  o.separator ''
  o.on("-r", "--rec-only", "Run only the reconstruction, then the analysis benchmark") do |a|
    opt.run_rec_down = true
    required_set = true
  end
  o.separator ''
  o.on("-b", "--ana-only", "Run only the analysis benchmark") do |a|
    opt.run_ana_down = true
    required_set = true
  end
  o.separator ''
  o.on("-s", "--sim-mode [SIMULATION_MODE]", "Run the simulation, reconstruction, and analysis",
       "[SIMULATION_MODE] must be one of:") do |a|
    unless avail_sim_modes.include? a
      $stderr.puts "ERROR: unknown simulation mode '#{a}'"
      exit 1
    end
    opt.sim_mode = a
    required_set = true
  end
  avail_sim_modes.each{ |it| o.separator ' '*40+it }
  o.separator ''
  o.separator 'optional options:'.upcase
  o.separator ''
  o.on("--sim-file [FILE]", "simulation file name",     "default = #{opt.sim_file}") { |a| opt.sim_file=a }
  o.on("--rec-file [FILE]", "reconstruction file name", "default = #{opt.rec_file}") { |a| opt.rec_file=a }
  o.on("--ana-file [FILE]", "analysis file name",       "default = #{opt.ana_file}") { |a| opt.ana_file=a }
  o.separator ''
  o.on("--[no-]run-sim", "simulation on/off",         "default = #{opt.run_sim ? 'on' : 'off'}") { |a| opt.run_sim=a }
  o.on("--[no-]run-rec", "reconstruction on/off",     "default = #{opt.run_rec ? 'on' : 'off'}") { |a| opt.run_rec=a }
  o.on("--[no-]run-ana", "analysis benchmark on/off", "default = #{opt.run_ana ? 'on' : 'off'}") { |a| opt.run_ana=a }
  o.separator ''
  o.on("-n", "--num-events [NUM_EVENTS]", Integer, "Number of events", "default = #{opt.num_events}") { |a| opt.num_events=a }
  o.on("-a", "--algos [ALGORITHMS]...", Array, "List of analysis algorithms to run; default = all",
       "delimit by commas, no spaces",
       "for more info, run:  #{opt.benchmark_exe}") { |a| opt.algos=a }
  o.on("-p", "--plots [PLOTS]...",      Array, "List of plots to draw", 
       "delimit by commas, no spaces",
       "default = all") { |a| opt.plots=a }
  o.on("-x", "--benchmark-exe [EXECUTABLE]", "benchmark executable",
       "default = #{opt.benchmark_exe} (from $PATH)") { |a| opt.benchmark_exe=a}
  o.separator ''
  o.on("-v", "--verbose", "Increase verbosity (-vv for more verbose)") { |a| opt.verbosity+=1 }
  o.on("-d", "--dry-run", "Dry run (just print commands)") { |a| opt.dry_run=true }
  o.on("--ci", "output plots to ./results, for CI artifact collection") { |a| opt.using_ci=true }
  o.separator ''
  o.on_tail("-h", "--help", "Show this message") do
    puts o
    exit 2
  end
end.parse!(ARGV.length>0 ? ARGV : ['--help'])

# print options
puts 'settings: {'.upcase
opt.each_pair { |k,v| puts "#{k.to_s.rjust(20)} => #{v}" }
puts '}'

# check for required options
unless required_set
  $stderr.puts "ERROR: required options have not been set"
  $stderr.puts "run '#{$0} --help' for guidance"
  exit 1
end

# figure out which steps to run
run_step = { :sim=>false, :rec=>false, :ana=>false, }
if opt.run_ana_down
  run_step[:ana] = true
elsif opt.run_rec_down
  run_step[:rec] = true
  run_step[:ana] = true
else
  run_step[:sim] = opt.run_sim
  run_step[:rec] = opt.run_rec
  run_step[:ana] = opt.run_ana
end
puts "steps to run: #{run_step}"

# get compact file
if ENV['DETECTOR_PATH'].nil? or ENV['DETECTOR_CONFIG'].nil?
  $stderr.puts "ERROR: unknown DETECTOR_PATH or DETECTOR_CONFIG"
  exit 1
end
compact_file = "#{ENV['DETECTOR_PATH']}/#{ENV['DETECTOR_CONFIG']}.xml"


# simulation command generators
# ---------------------------------------------------
def theta2xyz(theta)
  [ Math.sin(theta), 0.0, Math.cos(theta) ]
end

def eta2theta(eta)
  2 * Math.atan(Math.exp(-eta))
end


# fixed angle particle gun
simulate_fixed_angle = Proc.new do |theta, energy, particle|
  [
    "npsim",
    "--runType batch",
    "--compactFile #{compact_file}",
    "--outputFile #{opt.sim_file}",
    "--part.userParticleHandler=''",  # allow opticalphotons in MC particles
    "--enableGun",
    "--numberOfEvents #{opt.num_events}",
    "--gun.particle #{particle}",
    "--gun.energy #{energy}*GeV",
    "--gun.direction '(#{theta2xyz(theta).join ", "})'",
  ]
end


# define simulation command
# ---------------------------------------------------
sim_cmd = Array.new
case opt.sim_mode
when /^fixedEta/
  key       = opt.sim_mode.sub('fixedEta','').downcase.to_sym
  fixed_eta = EtaTestValues[key]
  if fixed_eta.nil?
    $stderr.puts "ERROR: EtaTestValues[#{key}] is not defined"
    exit 1
  end
  fixed_theta = eta2theta fixed_eta
  puts """Simulating fixed-eta events:
  - eta      = #{fixed_eta}
  - theta    = #{180.0 * fixed_theta / Math::PI} degrees
  - energy   = #{IdealEnergy} GeV
  - particle = #{IdealParticle}
  """
  sim_cmd = simulate_fixed_angle.call fixed_theta, IdealEnergy, IdealParticle
else
  exit 1 if run_step[:sim]
end


# define reconstruction command
# ---------------------------------------------------
output_collections = [
  "DRICHHits",
  "MCParticles",
  "DRICHRawHits",
  "DRICHRawHitsAssociations",
  "DRICHAerogelTracks",
  "DRICHGasTracks",
  "DRICHAerogelIrtCherenkovParticleID",
  "DRICHGasIrtCherenkovParticleID",
  "DRICHMergedIrtCherenkovParticleID",
  "ReconstructedChargedParticleAssociationsWithDRICHPID",
]
recon_cmd = [
  'eicrecon',
  "-Ppodio:output_collections=\"#{output_collections.join ','}\"",
  '-Pjana:nevents="0"',
  '-Pjana:debug_plugin_loading="1"',
  '-Pacts:MaterialMap="calibrations/materials-map.cbor"',
  "-Ppodio:output_file=\"#{opt.rec_file}\"",
  '-PDRICH:DRICHIrtCherenkovParticleID:cheatPhotonVertex=true',  # allow knowledge of true photons, for accurate residual determination
  opt.sim_file,
]


# define analysis benchmark command
# ---------------------------------------------------
analysis_cmd = [
  opt.benchmark_exe,
  "-i #{opt.rec_file}",
  "-o #{opt.ana_file}",
]
analysis_cmd.append "-a #{opt.algos.join ' '}" if opt.algos.size > 0
analysis_cmd.append '-' + 'v'*opt.verbosity if opt.verbosity > 0

# define analysis draw command
# ---------------------------------------------------
draw_cmd = [
  "#{__dir__}/draw_benchmark.py",
  "-i #{opt.ana_file}",
  "-o #{opt.using_ci ? "results/#{opt.sim_mode}" : opt.ana_file.gsub(/edm4hep.root$/,"plots")}"
]

# execute commands
# ---------------------------------------------------

# proc: execute a command; raise runtime exception if failure (exit nonzero)
exe = Proc.new do |cmd_args, name, step|
  if run_step[step]
    case step
    when :sim
      FileUtils.mkdir_p File.dirname(opt.sim_file)
    when :rec
      FileUtils.mkdir_p File.dirname(opt.rec_file)
    when :ana
      FileUtils.mkdir_p File.dirname(opt.ana_file)
    end
    cmd = cmd_args.join ' '
    puts "benchmark #{name} command:".upcase
    cmd_args.each_with_index do |arg,i|
      line = i==0 ? '' : '  '
      line += arg
      line += ' \\' unless i+1==cmd_args.size
      puts line
    end
    unless opt.dry_run
      puts '-'*50
      puts "#{name} execution:".upcase
      system cmd or raise "benchmark #{name} failed!".upcase
    end
    puts '-'*50
  end
end
puts '-'*50

# execute the commands
exe.call sim_cmd,      'simulation',     :sim
exe.call recon_cmd,    'reconstruction', :rec
exe.call analysis_cmd, 'analysis',       :ana
exe.call draw_cmd,     'draw',           :ana
