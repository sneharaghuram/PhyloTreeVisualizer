// components/ui/button.tsx
import * as React from "react";

// Use built-in typing so all native <button> props are supported
export type ButtonProps = React.ButtonHTMLAttributes<HTMLButtonElement>;

const Button = React.forwardRef<HTMLButtonElement, ButtonProps>(
  ({ children, className = "", type = "button", disabled = false, ...props }, ref) => {
    return (
      <button
        ref={ref}
        type={type}
        disabled={disabled}
        className={`btn inline-flex items-center justify-center rounded-md border px-3 py-2 text-sm font-medium ${className}`}
        {...props}
      >
        {children}
      </button>
    );
  }
);

Button.displayName = "Button";
export default Button;
