require "test/unit"

dir = File.dirname(__FILE__)

Test::Unit::AutoRunner.run(true, File.join(dir, "lin"))

Test::Unit::AutoRunner.run(true, File.join(dir, "eig"))
