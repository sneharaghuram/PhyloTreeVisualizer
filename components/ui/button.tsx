import React from "react";

// Define ButtonProps interface for typing the props
interface ButtonProps {
  children: React.ReactNode;  // Button text or elements passed as children
  onClick: () => void;        // Function to handle button click
  className?: string;         // Optional custom className
  type?: "button" | "submit" | "reset"; // Optional, default is "button"
  disabled?: boolean;         // Optional, default is false
}

const Button: React.FC<ButtonProps> = ({
  children,
  onClick,
  className = "",
  type = "button",
  disabled = false,
}) => {
  return (
    <button
      type={type}
      onClick={onClick}
      className={`btn ${className}`}
      disabled={disabled}
    >
      {children} {/* Button text or elements */}
    </button>
  );
};

export default Button;
