#ifndef SDMMSTOOLARGEEXCEPTION
#define SDMMSTOOLARGEEXCEPTION

namespace casac 
{
    //
    // A class to describe SDMMSTooLarge Exception.
    //
    class SDMMSTooLargeException {  
      public:
        /**
         * An empty contructor.
         */
        SDMMSTooLargeException();
  
        /**
         * A constructor with a message associated with the exception.
         * @param m a string containing the message.
         */
        SDMMSTooLargeException(string m);
  
        /**
         * The destructor.
         */
        virtual ~SDMMSTooLargeException();
  
        /**
         * Returns the message associated to this exception.
         * @return a string.
         */
        string getMessage() const;

      protected:
        string message;
    };

    inline SDMMSTooLargeException::SDMMSTooLargeException() : message ("SDMMSTooLargeException") {}
    inline SDMMSTooLargeException::SDMMSTooLargeException(string m) : message(m) {}
    inline SDMMSTooLargeException::~SDMMSTooLargeException() {}
    inline string SDMMSTooLargeException::getMessage() const {
        return "SDMMSTooLargeException : " + message;
    }
}

#endif
