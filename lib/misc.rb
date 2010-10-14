# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

module Misc
  def log(msg)
    $stderr.puts "#{Time.new}: #{msg}"
  end

  def usage(msg=nil, ec=1)
    printf "\nUPS!: #{msg}\n\n" if msg
    puts DATA.read
    exit(ec)
  end
end

