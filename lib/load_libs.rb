# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

bin_dir  = File.dirname($0)
main_dir = File.dirname(bin_dir)

lib_dir = File.join(main_dir, "lib")
Dir[File.join(lib_dir, "*.rb")].each do |file| 
  require File.basename(file)
end

# Make sure we use ruby >= 1.9
min_release  = "1.9.0"
ruby_release = "#{RUBY_VERSION}"

unless ruby_release.split('.')[0].to_i >= 1 &&
       ruby_release.split('.')[1].to_i >= 9
  $stderr.puts "This software requires a ruby version >= #{min_release}"
  $stderr.puts "This tool was run using: #{ruby_release}"
  $stderr.puts "Do not panic, it is just a matter of upgrading your ruby version"
  exit 1
end
